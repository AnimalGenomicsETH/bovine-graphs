#!/usr/bin/env python

from flask import Flask, render_template, request, url_for, redirect
import graph_viz_ht
import os
from pathlib import Path
import glob
import argparse

app = Flask(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--graphtype", help="Graph to visualize")
    return parser.parse_args()


@app.route("/")
def frontpage():
    return render_template("layout.html")


@app.route("/viz", methods=["GET", "POST"])
def viz():
    # Don't litter, remove graph constructed previously
    if glob.glob("static/*_*"):
        for file in glob.glob("static/*_*"):
            os.remove(file)
    if request.method == "POST":
        start_node = request.form["start_node"]
        stop_node = request.form["stop_node"]
        #nodeinf[nodeid] = [int(rrank), int(nodelen), chromo, int(pos)]
        chromo = nodeinf[start_node][2]
        locat = nodeinf[start_node][3]
        graph_viz_ht.visualize_graph(graphtype, graphcomp, graph, nodeinf, start=start_node, stop=stop_node, outf="svg")
        file_path = url_for("static", filename=f"graph5_{start_node}_{stop_node}.svg")
        return render_template("viz_show.html",
                               start=start_node,
                               stop=stop_node,
                               chromo=chromo,
                               locat=locat,
                               file_path=file_path)
    # If refresh just go into the frontpage
    if request.method == "GET":
        return redirect(url_for('frontpage'))


if __name__ == "__main__":
    #graphtype = "graph5"
    #graphcomp = ["UCD", "Angus", "Highland", "Brahman", "Yak"]
    args = parse_args()
    graphtype = args.graphtype
    # get the component of graphtype
    with open("../config/graph_comp.tsv") as infile:
        for line in infile:
            line_comp = line.strip().split()
            if line_comp[0] == graphtype:
                graphcomp = line_comp[1].strip().split(",")
                break
            else:
                raise ValueError("Specified graph cannot be found")
    # get the working directory of the results
    dirname = Path(__file__).absolute().parents[1]
    with open("../config/config.yaml") as infile:
        for line in infile:
            line_comp = line.strip().split()
            if "workdir" in line:
                outname = line_comp[1]
    basegraph = str(dirname.joinpath(outname).joinpath("graph"))
    graph = graph_viz_ht.generate_edges(basegraph, graphtype)
    nodeinf = graph_viz_ht.graph_info(basegraph, graphtype)

    # app.run(debug=True)
    app.run()
