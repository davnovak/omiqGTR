# omiqgtR

This package offers **<a href="https://www.r-project.org">R</a>-<a href="https://www.omiq.ai">OMIQ</a> interoperability for gating**.
This means that a gating tree from OMIQ can be downloaded and applied to a flow/mass cytometry expression matrix using a single `gate` function in R.
This is different from exporting cell-wise flags, and allows you to gate previously unseen data.

*Work in progress, see the Limitations section below.*

### Quick start

**Step 1**. Download your gating hierarchy: open the **Gating** task in your OMIQ workflow and press Ctrl+Shift+D.

**Step 2**. Install `omiqgtR` in your R console: `pak::pkg_install("davnovak/omiqqtR")`.

**Step 3**. Import the gating tree and plot it in R: `gt <- omiqgtR::parse_omiqgt("path_to_gatingfile.omiqgt"); plot(gt)`.

**Step 4**. Apply the gating to an FCS file in R: `gating_matrix <- omiqgtR::gate(gt, "path_to_fcsfile.fcs")`.
Each row of the resulting matrix will correspond to a cell, and each (named) column to a gate.

(**Step 4.**) If your gating hierarchy contains per-file adjustments, be sure to specify the FCS file's *OmiqID* using the `omiq_id` argument of `omiqgtR::gate()`.
This can be found in the OMIQ Dataset metadata table (**Dataset** -> **File Metadata**).

### Troubleshooting

Check function help in R: `?parse_omiqgt`, `?print.GatingTree`, `?plot.GatingTree`, and `?gate`.

### Limitations

* Only rectangular and polygonal gates are currently supported.
* I have yet to implement interoperability with `flowWorkspace`, `CytoML`, or maybe even *FlowJo*.
If you want this, post a feature request in the *Issue* tab of this GitHub repository.
