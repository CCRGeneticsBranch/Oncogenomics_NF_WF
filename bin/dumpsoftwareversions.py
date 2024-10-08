#!/usr/bin/env python

import yaml
import platform
from textwrap import dedent
import sys
from datetime import datetime


def _make_versions_html(versions):
    html = [
        dedent(
            """\\
            <style>
            #nf-core-versions tbody:nth-child(even) {
                background-color: #f2f2f2;
            }
            </style>
            <table class="table" style="width:100%" id="nf-core-versions">
                <thead>
                    <tr>
                        <th> Process Name </th>
                        <th> Software </th>
                        <th> Version  </th>
                    </tr>
                </thead>
            """
        )
    ]
    for process, tmp_versions in sorted(versions.items()):
        html.append("<tbody>")
        for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
            html.append(
                dedent(
                    f"""\\
                    <tr>
                        <td><samp>{process if (i == 0) else ''}</samp></td>
                        <td><samp>{tool}</samp></td>
                        <td><samp>{version}</samp></td>
                    </tr>
                    """
                )
            )
        html.append("</tbody>")
    html.append("</table>")
    return "\\n".join(html)


task_process = sys.argv[3]
versions_this_module = {}
versions_this_module[task_process] = {
    "python": platform.python_version(),
    "yaml": yaml.__version__,
}

with open(sys.argv[1]) as f:
    versions_by_process = yaml.load(f, Loader=yaml.BaseLoader) | versions_this_module

# aggregate versions by the module name (derived from fully-qualified process name)
versions_by_module = {}
versions_by_module["Pipeline"] = {
    "version": sys.argv[5],
}
for process, process_versions in versions_by_process.items():
    module = process.split(":")[-1]
    try:
        assert versions_by_module[module] == process_versions, (
            "We assume that software versions are the same between all modules. "
            "If you see this error-message it means you discovered an edge-case "
            "and should open an issue in nf-core/tools. "
        )
    except KeyError:
        versions_by_module[module] = process_versions

versions_by_module["Nextflow"] = {
    "Nextflow": sys.argv[2],
}

timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
Patient = sys.argv[4]
output_file = f"{Patient}.config.{timestamp}.txt"

with open(output_file, "w") as f:
    yaml.dump(versions_by_module, f, default_flow_style=False)
