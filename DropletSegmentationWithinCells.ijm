/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU Affero General Public License version 3 as published by the Free Software Foundation:
 http://www.gnu.org/licenses/agpl-3.0.txt
*/

// input parameters
#@ File (label = "Input directory", style = "directory") input
#@ Integer (label = "Droplet channel", value = 1, style = "spinner") droplet_channel
#@ Integer (label = "Nucleus channel", value = -1, style = "spinner") nucleus_channel
#@ Integer (label = "Cell size for Cellpose", value = 100, style = "spinner") cell_diameter
#@ Float (label = "Threshold for Cellpose", value = 0.0, style="format:#.##") cell_prob_threshold
#@ Integer (label = "Minimum cell size", value = 1000, style = "spinner") min_cell_size
#@ Float (label = "LoG radius for droplets", value = 2.0, style="format:#.##") LoG_radius
#@ Float (label = "LoG quality for droplets", value = 50.0, style="format:#.##") LoG_quality
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix

// call to the main function "processFolder"
processFolder(input);

// function to scan folders to find files with correct suffix
function processFolder(input) {
	///////////// initial cleaning /////////////////
	// close all images
	run("Close All");
	// reset ROI manager
	roiManager("Reset");
	// clear results
	run("Clear Results");

	///////////// apply pipeline to input images /////////////////
	// get the files in the input folder
	list = getFileList(input);
	list = Array.sort(list);
	// loop over the files
	for (i = 0; i < list.length; i++) {
		// if current file ends with the suffix given as input parameter, call function "processFile" to process it
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
	// save parameters
	// create results table
	Table.create("Results");
	setResult("Droplet channel", 0, droplet_channel);
	setResult("Nucleus channel", 0, nucleus_channel);
	setResult("Cell size for Cellpose", 0, cell_diameter);
	setResult("Threshold for Cellpose", 0, cell_prob_threshold);
	setResult("Minimum cell size", 0, min_cell_size);
	setResult("LoG radius for droplets", 0, LoG_radius);
	setResult("LoG quality for droplets", 0, LoG_quality);
	updateResults();
	// save results
	saveAs("Results", output + File.separator + "parameters.csv");	
	// close results table
	close("Results"); 
	// reset ROI manager
	roiManager("Reset");

}

function processFile(input, output, file) {

	// open image
	open(input + File.separator + file);
	// rename
	rename("Input");
	
	// segment cells with cellpose
	run("Cellpose Advanced", "diameter=" + cell_diameter + " cellproba_threshold=" + cell_prob_threshold + " flow_threshold=10.0 anisotropy=1.0 diam_threshold=12.0 model=cyto2 nuclei_channel=" + nucleus_channel + " cyto_channel=" + droplet_channel + " dimensionmode=2D stitch_threshold=0 omni=false cluster=false additional_flags=4");
	// remove cells on the border
	run("Remove Border Labels", "left right top bottom");
	// size filtering
	run("Label Size Filtering", "operation=Greater_Than size=" + min_cell_size + "");
	// connected components
	run("Connected Components Labeling", "connectivity=4 type=[16 bits]");
	// size filtering second round as isolated pixels can be considered
	run("Label Size Filtering", "operation=Greater_Than size=" + min_cell_size + "");
	// connected components
	run("Connected Components Labeling", "connectivity=4 type=[16 bits]");
	rename("Cells");
	
	// get image dimensions
	// select input image
	selectImage("Input");
	getDimensions(width, height, nChannels, nbSlices, frames);
	
	if (nChannels == 1) {
		// select input image
		selectImage("Input");
		// run LoG
		run("LoG Trackmate", "radius=" + LoG_radius + " quality=" + LoG_quality + "");
	}
	else {
		// select input image
		selectImage("Input");
		// split channels
		run("Split Channels");
		// select droplet channel
		selectImage("C" + droplet_channel + "-Input");
		// change image properties
		run("LoG Trackmate", "radius=" + LoG_radius + " quality=" + LoG_quality + "");
	}
	// connected components
	run("Connected Components Labeling", "connectivity=4 type=[16 bits]");
	rename("DropletCenters");
	// dilation
	// segmented droplets
	run("Morphological Filters", "operation=Dilation element=Disk radius=" + LoG_radius + "");
	rename("Droplets");
	
	// get number of cells
	selectImage("Cells");
	getStatistics(area, mean, min, nb_cells, std, histogram);
	// get number of droplets per cell
	nb_droplets_per_cell = newArray(nb_cells);
	// compute intensity in cell image at each droplet location
	run("Intensity Measurements 2D/3D", "input=Cells labels=DropletCenters mean");
	// rename measurements as Results
	Table.rename("Cells-intensity-measurements", "Results");
	// get cell id associated to each droplet -> #droplets per cell
	for (i = 0; i < nResults; i++) {
		cell_id = getResult("Mean", i);
		if (cell_id>0){
			nb_droplets_per_cell[cell_id-1] += 1;
		}
	}
	// close results table
	close("Results"); 

	// visual inspection of results
	if (nChannels > 1) {
		// merge channels	
		run("Merge Channels...", "c1=C" + droplet_channel + "-Input c2=C" + nucleus_channel + "-Input create");
		rename("Input");
	}
	// select droplet images and get ROIs
	selectImage("Droplets");
	run("Label image to ROIs");
	// overlay droplets on input image
	selectImage("Input");
	roiManager("Show All without labels");
	run("Flatten");
	rename("ResultImage");
	// select cells image and get ROIs
	selectImage("Cells");
	run("Label image to ROIs");
	// get cells area
	cells_area = newArray(nb_cells);
	for (i = 0; i < nb_cells; i++) {
		roiManager("select", i);
		getStatistics(cells_area[i], mean, min, max, std, histogram);
	}
	// overlay cells on result image with their labels
	selectImage("ResultImage");
	roiManager("Show All with labels");
	run("Flatten");
	// save for visual inspection
	saveAs("png", output + File.separator + file);
	
	// create result output with number of droplets per cells and cell area
	for (i = 0; i < nb_cells; i++) {
		setResult("Cell id", i, i+1);
		setResult("#droplets", i, nb_droplets_per_cell[i]);
		setResult("Cell area", i, cells_area[i]);	
	}
	updateResults();	
	// save results
	file_name = replace(file, ".tif", "_results.csv");
	saveAs("Results", output + File.separator + file_name);	

	///////////// clear everything /////////////////
	// close all images
	run("Close All");
	// reset ROI manager
	roiManager("Reset");
	// close results table
	close("Results"); 

}
