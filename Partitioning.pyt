#----------------------------------------------------------------------
# Name:        ibtool
# Purpose:     Toolset fuer die Siedlungsabgrenzung auf Grundlage von
#              Gebaeudegrundrissen
#
# Author:      OLiver Harig
#
# Created:     17.04.2014
# Copyright:   (c) Oliver Harig 2021
# Licence:     Alle Rechte vorbehalten
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
        self.label =  "ibtool Toolbox"
        self.alias  = "ibtool"
        self.description = "In diesem Toolset sind sowohl alle beoetigten Werkzeuge zu Erstellung " + \
                           "von Siedlungsgrenzen auf Grundlage von Gebaedegrundrissen als auch die Werkzeuge " + \
                           " zur Erstellung und von Abgrenzungen mit Hilfe von Trainingsdaten enthalten."

        # List of tool classes associated with this toolbox
        self.tools = [Partitionen]

class ibtool:

    import math
    arcpy.env.overwriteOutput = True



# ###################### SIEDLUNGSGRENZEN-TOOL ########################

    def siedgr(self, input_HU, cell_size, density_value):

        arcpy.AddMessage("Datenaufbereitung: Siedlungsgrenzen erzeugen")

        #SteuerVariablen
        radius = 2 * cell_size
        d_var = str(cell_size // 2)
        radius_del = '{} Meters'.format(d_var)


        #Datei-Aliase
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

        # laedt ArcGIS Spatial Analyst license
        arcpy.CheckOutExtension("Spatial")

        # Ueberschreiben = true
        arcpy.env.overwriteOutput = True

        arcpy.management.FeatureToPoint(Input_feature, Input_feature_point, "INSIDE")


        # erzeugt Raster mit Dichtewerten
        HU_dest = arcpy.sa.PointDensity(Input_feature_point, "NONE", cell_size, "Circle {} MAP".format(radius), "SQUARE_METERS");

        #out_raster.save(HU_dest)

        # Raster zu Punkt
        arcpy.RasterToPoint_conversion(HU_dest, HU_Raster, 'Value')

        # erzeugt Punktrastercluster
        arcpy.MakeFeatureLayer_management(HU_Raster, HU_Raster_Layer, '"grid_code" > {}'.format(density_value))
        arcpy.CopyFeatures_management(HU_Raster_Layer, HU_Raster_Feature)

        # erzeugt Thiessenpolygone
        arcpy.CreateThiessenPolygons_analysis(HU_Raster_Feature, Thiess_Poly, 'ONLY_FID')

        # Polygone zu Linien
        arcpy.FeatureToLine_management(Thiess_Poly, Thiess_Line, '#', 'ATTRIBUTES')

        # teilt die Thisssenlinien an den Schnittpunkten
        arcpy.SplitLine_management(Thiess_Line, Thiess_Split)

        # loescht alle Linien in der Naehe der Punktrastercluster
        arcpy.MakeFeatureLayer_management(Thiess_Split, Thiess_Split_Lay)
        arcpy.SelectLayerByLocation_management(Thiess_Split_Lay, 'WITHIN_A_DISTANCE', HU_Raster_Layer, radius_del, 'NEW_SELECTION')
        arcpy.DeleteFeatures_management(Thiess_Split_Lay)

        # erzeugt Siedlungsgrenzenpolygon
        arcpy.FeatureToPolygon_management(Thiess_Split_Lay, Poly_Grenz, '#', 'ATTRIBUTES', '#')

        # loescht die TempDateien
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
        self.label       = "Datensatzpartionierung"
        self.description = "Mit diesem Tool werden aus gegeben Gebaeudegrundrissen " + \
                           "Partitionierungspolygone fuer die schrittweise Abarbeitung des Algorithmus erstellt."

        self.canRunInBackground = False

    def getParameterInfo(self):

        param0 = arcpy.Parameter(
            displayName="Workspace",
            name="arcpy.env.workspace",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        param1 = arcpy.Parameter(
            displayName="Gebauedegrundrisse",
            name="imput_HU",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input")
        param1.filter.list = ["POLYGON"]

        param2 = arcpy.Parameter(
            displayName="Zellengroesse (Meter)",
            name="cell_size",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        param2.value = 150

        param3 = arcpy.Parameter(
            displayName="Dichteschwellwert",
            name="density_value",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param3.value = 0.00001


        param4 = arcpy.Parameter(
            displayName="Ausgabepolygon",
            name="out_put_poly",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")
        param4.filter.list = ["POLYGON"]

        param5 = arcpy.Parameter(
            displayName="Referenzabgrenzung",
            name="Muster_grenz",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input")
        param5.category = "Konfliktpruefung (optional)"

        params = [param0, param1, param2, param3, param4, param5]

        return params

    def execute(self, params, messages):

            #  Lokale Variablen
            arcpy.env.workspace = params[0].valueAsText
            input_HU = params[1].valueAsText
            cell_size = params[2].valueAsText
            density_value = params[3].valueAsText
            out_put_poly = params[4].valueAsText
            Muster_grenz = params[5].valueAsText


            #  Hauptprogramm
            out_siedgr = ibtool.siedgr(input_HU, int(cell_size), float(density_value))

            #  Konfliktpruefung
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


            # Ausgabe
            arcpy.CopyFeatures_management(out_siedgr, out_put_poly)

