import ROOT

run_name_list = ["r1159", "r0990", "r1010", "r0948"]

base_directory = "/media/olivia/Partition1/CATS"
file_pattern_by_run = {run_name: f"{base_directory}/*{run_name}*.root" for run_name in run_name_list}

tree_name = "AD"
draw_expression = "CATS1XVN"
cut_expression = "CATS1XV>300"
draw_options = ""

canvas_columns = 2
canvas = ROOT.TCanvas("canvas_runs", "CATS1XVN per run", 1400, 900)
number_of_rows = (len(run_name_list) + canvas_columns - 1) // canvas_columns
canvas.Divide(canvas_columns, number_of_rows)

chain_by_run = {}
histogram_by_run = {}

for pad_index, run_name in enumerate(run_name_list, start=1):
    chain_object = ROOT.TChain(tree_name)
    number_of_files_added = chain_object.Add(file_pattern_by_run[run_name])
    chain_by_run[run_name] = chain_object

    canvas.cd(pad_index)

    histogram_name = f"h_{draw_expression}_{run_name}"
    ROOT.gDirectory.Delete(histogram_name + ";*")

    chain_object.Draw(f"{draw_expression}>>{histogram_name}", cut_expression, draw_options)
    histogram_object = ROOT.gDirectory.Get(histogram_name)
    histogram_by_run[run_name] = histogram_object

    if histogram_object:
        histogram_object.SetTitle(f"{run_name}  (files: {number_of_files_added}, entries: {chain_object.GetEntries()})")
        histogram_object.Draw()
    else:
        empty_histogram = ROOT.TH1F(f"empty_{run_name}", f"{run_name} (no histogram)", 10, 0.0, 10.0)
        empty_histogram.Draw()

canvas.Update()
ROOT.gSystem.ProcessEvents()
print("DONE")
input("Press Enter to exit...")
