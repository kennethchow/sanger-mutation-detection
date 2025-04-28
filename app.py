import os
import subprocess
import tempfile
import yaml
from datetime import datetime
from markupsafe import Markup
from flask import (
    Flask, request, render_template, send_from_directory,
    flash, redirect, url_for
)

app = Flask(__name__)
app.secret_key = "replace-with-your-own-secret"

@app.template_filter("toyaml")
def to_yaml_filter(data):
    """
    Dump a Python object as YAML for use in a <textarea>.
    We mark it safe so the newlines show up properly.
    """
    dumped = yaml.safe_dump(data, sort_keys=False)
    # wrap in Markup so Jinja doesn't HTML-escape the newlines
    return Markup(dumped)

PIPELINE_SCRIPT = os.path.join(
    os.path.dirname(__file__),
    "sanger_blast_pipeline.py"
)

# keep track of last output_dir and csv file
LAST_OUTPUT_DIR = None
LAST_CSV         = None

# minimal config skeleton to prefill
DEFAULT_CFG = {
    "input_dir": "",
    "output_dir": "",
    "reference": {
        "fasta": "",
        "make_blast_db": True,
        "blast_db_name": "",
        "blast_db_dir": ""
    },
    "blast": {
        "task": "blastn-short",
        "outfmt": "5",
        "evalue": 1e-5
    },
    # example single mutation
    "mutations": {
        "L858R": {
            "type": "SNV",
            "pos": 177791,
            "wt": "T",
            "mut": "G"
        }
    }
}


@app.route("/", methods=["GET","POST"])
def index():
    global LAST_OUTPUT_DIR, LAST_CSV

    cfg = DEFAULT_CFG.copy()
    log = None
    csv_ready = False

    if request.method == "POST":
        # 1) collect fields
        cfg["input_dir"]  = request.form["input_dir"].strip()
        cfg["output_dir"] = request.form["output_dir"].strip()

        ref = cfg["reference"]
        ref["fasta"]        = request.form["reference_fasta"].strip()
        ref["make_blast_db"]= (request.form.get("make_blast_db")=="on")
        ref["blast_db_name"]= request.form["blast_db_name"].strip()
        ref["blast_db_dir"] = request.form["blast_db_dir"].strip()

        blast = cfg["blast"]
        blast["task"]   = request.form["blast_task"].strip()
        blast["outfmt"] = request.form["blast_outfmt"].strip()
        try:
            blast["evalue"] = float(request.form["blast_evalue"])
        except ValueError:
            flash("E-value must be a number", "danger")

        # 2) parse the YAML block for mutations
        muts_yaml = request.form["mutations_yaml"]
        try:
            muts = yaml.safe_load(muts_yaml)
            if not isinstance(muts, dict):
                raise ValueError("Top‐level must be a mapping")
            cfg["mutations"] = muts
        except Exception as e:
            flash(f"Mutations YAML error: {e}", "danger")

        # 3) write that config to a temp file
        os.makedirs(cfg["output_dir"], exist_ok=True)
        fd, cfg_path = tempfile.mkstemp(suffix=".yaml")
        os.close(fd)
        with open(cfg_path, "w") as f:
            yaml.safe_dump(cfg, f)

        # 4) run your pipeline
        cmd = ["python", PIPELINE_SCRIPT, "--config", cfg_path]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        log = proc.stdout + "\n" + proc.stderr

        if proc.returncode != 0:
            flash(f"Pipeline exited with code {proc.returncode}", "danger")

        # 5) locate the CSV and then rename it with timestamp + blast_db basename
        out_dir = cfg["output_dir"]
        # find the original CSV
        orig_csvs = [f for f in os.listdir(out_dir) if f.lower().endswith(".csv")]
        if orig_csvs:
            orig = orig_csvs[0]
            # build a new name
            tstamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            # use the blast_db_name from the config (or fallback)
            dbname = os.path.basename(cfg["reference"]["blast_db_name"]) or "egfr"
            new_name = f"{dbname}_{tstamp}.csv"
            # do the rename/move
            os.replace(
                os.path.join(out_dir, orig),
                os.path.join(out_dir, new_name)
            )
            LAST_OUTPUT_DIR = out_dir
            LAST_CSV         = new_name
            csv_ready        = True
        else:
            flash("No CSV found in output directory!", "warning")

        # keep the form filled
        # (cfg already updated)

    return render_template(
        "index.html",
        cfg=cfg,
        log=log,
        csv_ready=csv_ready,
        csv_filename=LAST_CSV
    )


@app.route("/download/<filename>")
def download(filename):
    if not LAST_OUTPUT_DIR or not filename:
        return "Nothing to download", 404
    return send_from_directory(LAST_OUTPUT_DIR, filename, as_attachment=True)


if __name__ == "__main__":
    app.run(debug=True)