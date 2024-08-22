# This program is free software; you can redistribute it and/or modify it under the terms of the GNU Affero General Public License version 3 as published by the Free Software Foundation:
# http://www.gnu.org/licenses/agpl-3.0.txt

#@ Float radius
#@ Float quality
import sys
 
from ij import IJ
from ij import WindowManager
from jarray import zeros
from ij.process import FloatProcessor
from ij import IJ, ImagePlus, ImageStack

from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
 
# We have to do the following to avoid errors with UTF8 chars generated in 
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')
 
# get current image
imp = IJ.getImage() 
 
#----------------------------
# Create the model object now
#----------------------------
 
# Some of the parameters we configure below need to have
# a reference to the model at creation. So we create an
# empty model now.

model = Model()
 
# Send all messages to ImageJ log window.
model.setLogger(Logger.IJ_LOGGER)
 
 
 
#------------------------
# Prepare settings object
#------------------------
 
settings = Settings(imp)
 
# Configure detector - We use the Strings for the keys
settings.detectorFactory = LogDetectorFactory()
settings.detectorSettings = {
    'DO_SUBPIXEL_LOCALIZATION' : True,
    'RADIUS' : radius,
    'TARGET_CHANNEL' : 1,
    'THRESHOLD' : 0.,
    'DO_MEDIAN_FILTERING' : False,
}  
 
# Configure spot filters - Classical filter on quality
filter1 = FeatureFilter('QUALITY', quality, True)
settings.addSpotFilter(filter1)
 
# Configure tracker - We want to allow merges and fusions
settings.trackerFactory = SparseLAPTrackerFactory()
settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = True
settings.trackerSettings['ALLOW_TRACK_MERGING'] = True
 
# Add ALL the feature analyzers known to TrackMate. They will 
# yield numerical features for the results, such as speed, mean intensity etc.
settings.addAllAnalyzers()
 
# Configure track filters - We want to get rid of the two immobile spots at
# the bottom right of the image. Track displacement must be above 10 pixels.
 
filter2 = FeatureFilter('TRACK_DISPLACEMENT', 10, True)
settings.addTrackFilter(filter2)
 
 
#-------------------
# Instantiate plugin
#-------------------
 
trackmate = TrackMate(model, settings)
 
#--------
# Process
#--------
 
ok = trackmate.checkInput()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))
 
ok = trackmate.process()
if not ok:
	IJ.log( 'No spots were detected' ) 
    
else:
 
	#----------------
	# Display results
	#----------------

	# A selection.
	selectionModel = SelectionModel( model )
 
	# Read the default display settings.
	ds = DisplaySettingsIO.readUserDefault()
 
	displayer =  HyperStackDisplayer( model, selectionModel, imp, ds )
	displayer.render()
	displayer.refresh()
 
	# Echo results with the logger we set at start:
	model.getLogger().log( str( model ) )

	# Get image dimensions
	width = imp.getWidth()
	height = imp.getHeight()
	depth = imp.getNSlices()

	# create stack
	stack = ImageStack(width, height)

	# extract spots
	spots = model.getSpots()

	# loop over z
	for i in range(depth):
		# create current slice
		pixels = zeros(width * height, 'f')
		# loop over spots, creation of a point when spot is in the current slice
		for spot in spots.iterator(True):
			sid = spot.ID()
			z=spot.getFeature('POSITION_Z')
			if int(z)==(i):
				x=spot.getFeature('POSITION_X')
				y=spot.getFeature('POSITION_Y')
				pixels[int(y)*width + int(x)] = int(sid)
		# current slice creation
		fp = FloatProcessor(width, height, pixels)
		# add slice to stack
		stack.addSlice(str(i), fp)

	# droplet center creation
	imp_centers = ImagePlus("LoGOutputCenters", stack)  
	imp_centers.show()
