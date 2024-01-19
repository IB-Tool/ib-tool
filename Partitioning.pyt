#----------------------------------------------------------------------
# Name:        IB-Tool2
# Purpose:     Toolset for the delineation of settlements on the basis
#              building footprints
#
# Author:      Oliver Harig
#
# Created:     10.05.2021
# Cite:        Oliver Harig (2021). Toolset for the delineation of settlements on the basis building footprints, road network and land use data (v1.0) https://doi.org/10.26084/IOERFDZ-SOFT-001    
# Licence:     MIT License
               Copyright 2021 Oliver Harig
#
#              Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
#              associated documentation files (the "Software"), to deal in the Software without restriction, including
#              without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#              copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the
#              following conditions:
#
#              The above copyright notice, and this permission notice must be included in all copies or substantial portions of the Software.
#
#              THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
#              LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#              IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#              WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#              SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#----------------------------------------------------------------------

import arcpy
import time
from arcpy import env
from numpy import*
import math
import random
from arcpy.sa import*

arcpy.env.referenceScale = "10000"
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

# ############################### PARTITIONING-TOOL ##########################

class Toolbox(object):
    def __init__(self):
        self.label =  "ib-tool2_partitioning "
        self.alias  = "partitioning "
        self.description = "ArcGIS toolbox for partitioning of geodata based on building footprints"

        # List of tool classes associated with this toolbox
        self.tools = [Partitionen]

class ibtool:

    arcpy.env.overwriteOutput = True

    def siedgr(self, input_HU, cell_size, density_value):

        arcpy.AddMessage("Start partitioning")

        #  Variables
        radius = 2 * cell_size
        d_var = str(cell_size // 2)
        radius_del = '{} Meters'.format(d_var)


        #  File-Aliases
        Input_feature = input_HU
        Input_feature_point= 'Input_feature_point'
        HU_Raster = 'HU_Raster'
        HU_Raster_Layer = 'HU_Raster_Layer'
        HU_Raster_Feature = 'HU_Raster_Feature'
        Thiess_Poly = 'Thiess_Poly'
        Thiess_Line = 'Thiess_Line'
        Thiess_Split = 'Thiess_Split'
        Thiess_Split_Lay = 'Thiess_Split_Lay'
        Poly_Grenz = 'Poly_Grenz_{}_{}'.format(cell_size, radius)

        #  loading ArcGIS Spatial Analyst license
        arcpy.CheckOutExtension("Spatial")

        arcpy.env.overwriteOutput = True
        arcpy.management.FeatureToPoint(Input_feature, Input_feature_point, "INSIDE")

        #  Creates rasters with density values
        HU_dest = arcpy.sa.PointDensity(Input_feature_point, "NONE", cell_size, "Circle {} MAP".format(radius), "SQUARE_METERS")

        #  Raster to point
        arcpy.RasterToPoint_conversion(HU_dest, HU_Raster, 'Value')

        #  Generates point grid clusters
        arcpy.MakeFeatureLayer_management(HU_Raster, HU_Raster_Layer, '"grid_code" > {}'.format(density_value))
        arcpy.CopyFeatures_management(HU_Raster_Layer, HU_Raster_Feature)

        #  generates Thiessen polygons
        arcpy.CreateThiessenPolygons_analysis(HU_Raster_Feature, Thiess_Poly, 'ONLY_FID')

        #  Polygons to lines
        arcpy.FeatureToLine_management(Thiess_Poly, Thiess_Line, '#', 'ATTRIBUTES')

        # divides the Thisssen lines at the intersections
        arcpy.SplitLine_management(Thiess_Line, Thiess_Split)

        # Deletes all lines in the proximity of the point grid clusters
        arcpy.MakeFeatureLayer_management(Thiess_Split, Thiess_Split_Lay)
        arcpy.SelectLayerByLocation_management(Thiess_Split_Lay, 'WITHIN_A_DISTANCE', HU_Raster_Layer, radius_del, 'NEW_SELECTION')
        arcpy.DeleteFeatures_management(Thiess_Split_Lay)

        # creates settlement boundary polygon
        arcpy.FeatureToPolygon_management(Thiess_Split_Lay, Poly_Grenz, '#', 'ATTRIBUTES', '#')

        # Delete the temporary files
        del_list = [HU_dest, HU_Raster, HU_Raster_Feature, HU_Raster_Layer, Thiess_Line, Thiess_Poly, Thiess_Split, Thiess_Split_Lay]
        for list in del_list:
            arcpy.Delete_management(list)

        #  Add Name-Field
        if len(arcpy.ListFields(Poly_Grenz, "NAME")) > 0:
            arcpy.DeleteField_management(Poly_Grenz, "NAME")
        arcpy.AddField_management(Poly_Grenz, "NAME", "TEXT")
        arcpy.management.CalculateField(Poly_Grenz, "NAME", "'PART_'+ str(!OBJECTID!)", "PYTHON_9.3", None)

        return Poly_Grenz



ibtool = ibtool()


class Partitionen(object):

    def __init__(self):
        self.label       = "Data set partitioning"
        self.description = "With this tool, partitioning polygons are created from given building footprints for the step-by-step processing of the algorithm."

        self.canRunInBackground = False

    def getParameterInfo(self):

        param0 = arcpy.Parameter(
            displayName="Workspace",
            name="arcpy.env.workspace",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        param1 = arcpy.Parameter(
            displayName="Building footprints",
            name="imput_HU",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input")
        param1.filter.list = ["POLYGON"]

        param2 = arcpy.Parameter(
            displayName="Cell size (metres)",
            name="cell_size",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        param2.value = 150

        param3 = arcpy.Parameter(
            displayName="Density threshold",
            name="density_value",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param3.value = 0.00001


        param4 = arcpy.Parameter(
            displayName="Output polygon",
            name="out_put_poly",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")
        param4.filter.list = ["POLYGON"]

        param5 = arcpy.Parameter(
            displayName="Expert delineation",
            name="Muster_grenz",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input")
        param5.category = "Conflict check (optional)"

        params = [param0, param1, param2, param3, param4, param5]

        return params

    def execute(self, params, messages):

            #  Local variables
            arcpy.env.workspace = params[0].valueAsText
            input_HU = params[1].valueAsText
            cell_size = params[2].valueAsText
            density_value = params[3].valueAsText
            out_put_poly = params[4].valueAsText
            Muster_grenz = params[5].valueAsText

            #  Main programme
            out_siedgr = ibtool.siedgr(input_HU, int(cell_size), float(density_value))

            #  Conflict check
            if str(Muster_grenz) != 'None':
                arcpy.FeatureToLine_management(out_siedgr, 'out_siedgr_line')
                arcpy.MakeFeatureLayer_management('out_siedgr_line','out_siedgr_line_FL')
                select_DELI = arcpy.SelectLayerByLocation_management('out_siedgr_line_FL', 'INTERSECT', Muster_grenz, '#', 'NEW_SELECTION')
                arcpy.CopyFeatures_management(select_DELI, 'DELI_sel')
                result = arcpy.GetCount_management('DELI_sel')
                ANZ = int(result.getOutput(0))
                arcpy.AddMessage("Es wurden " + str(ANZ) + " Konflikte mit der Referenzabgrenzung festgestelt.")
                arcpy.DeleteFeatures_management('DELI_sel')
                arcpy.DeleteFeatures_management('out_siedgr_line')

            #  Output
            arcpy.CopyFeatures_management(out_siedgr, out_put_poly)

