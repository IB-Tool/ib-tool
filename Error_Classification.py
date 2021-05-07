
# ------ IMPORTS ------

import arcpy
import os
import time

# ------ GLOBAL SETTINGS ------

# working folder set to current folder of py-file
Workspace = os.getcwd()
arcpy.env.overwriteOutput = True

# deactivate ESRI logging for disk space saving
arcpy.SetLogHistory(False)

# create folder for temporary files
if arcpy.Exists("Tmp.gdb"):
    arcpy.Delete_management("Tmp.gdb")
arcpy.CreateFileGDB_management(Workspace, "Tmp.gdb")
Geodatabase = Workspace + os.path.sep + "Tmp.gdb"
timedic = {}


# ------ INPUT DATA ------

InputHU = Workspace + os.path.sep + "A_HU_P.shp" #  Building footprints (polygon, shape)
IBS = Workspace + os.path.sep + "A_IBS_P.shp" #  Expert delieation (polygon, shape)
UGB = Workspace + os.path.sep + "A_UGB_P5.shp"  #  Caculated boundary (polygon, shape)
Nutzungen = Workspace + os.path.sep + "A_Nutzungen_P.shp" #  land use geometry (polygon, shape)
sr = arcpy.SpatialReference(25833)


# set thresholds

GOT = 15  # Upper threshold
LBC = 3   # Lower threshold

# set ESRI environment variables
arcpy.env.workspace = Workspace
arcpy.env.outputCoordinateSystem = sr
arcpy.env.referenceScale = "10000"
arcpy.env.outputCoordinateSystem = sr
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
os.chdir(Workspace)


# ------ GLOBAL FUNCTIONS ------

def tmp(FileName):
    # write file in temporary folder and return FileName
    return Workspace + os.path.sep + "Tmp" + os.path.sep + FileName


def mem(FileName):
    # write file in RAM-memory and return FileName
    return r'in_memory' + os.path.sep + FileName


def gdb(FileName):
    # write file in geodatabase and return FileName
    return Geodatabase + os.path.sep + FileName


def DelName(DelList):
    # delete files given as list of files
    arcpy.env.overwriteOutput = True
    if DelList != None:
        for j in DelList:
            try:
                arcpy.Delete_management(j)
            except Exception:
                pass


def DelField(Filename, Fields):
    # delete field of given file
    try:
        for i in Fields:
            arcpy.DeleteField_management(Filename, i)
    except Exception:
        pass


def TimePrint3(T3, text=""):
    # logging the time for the section for debugging
    T2 = time.time()
    Ta = (T2 - T3)
    Tdiff = 0

    if "{}_Ts".format(text) in timedic:
        timedic["{}_Tj".format(text)] = Ta
        Tdiff = timedic["{}_Tj".format(text)] - timedic["{}_Ts".format(text)]
    else:
        timedic["{}_Ts".format(text)] = Ta
    printtext = text + ": " + str(round(Ta, 1)) + ", " + str(round(Tdiff, 2))
    print (printtext)
    return T2


def Shp_Area(FileName, Fieldname='Shape_Area'):
    # adds shape area field to file
    if len(arcpy.ListFields(FileName, Fieldname)) == 0:
        arcpy.AddField_management(FileName, '{}'.format(Fieldname), "DOUBLE")
    arcpy.management.CalculateField(FileName, '{}'.format(Fieldname), "!shape.geodesicArea@SQUAREMETERS!", "PYTHON_9.3",
                                    None)


def Shp_Length(FileName):
    # adds shape length field to file
    if len(arcpy.ListFields(FileName, "Shape_Len")) == 0:
        arcpy.management.AddField(FileName, "Shape_Len", "DOUBLE")
    arcpy.management.CalculateField(FileName, "Shape_Len", "!shape.geodesicLength@METERS!", "PYTHON_9.3", None)


def CheckFileType(inputfile, datatype, geometrytype):
    """
    :param inputfile: File to check
    :param datatype (str): ShapeFile, FeatureClass, FeatureLayer
    :param geometrytype (str):  Polygon, Polyline, Point, Multipoint, MultiPatch

   -  checks whether the file has the required data type and geometry type
    """

    desc = arcpy.Describe(inputfile)
    if datatype != desc.dataType:  # ShapeFile FeatureClass FeatureLayer
        print ("Data type of {} should be {}, but is {}".format(inputfile, datatype, desc.dataType))
        return inputfile
    if geometrytype != desc.shapeType:
        print ("Geometry type of {} should be {}, but is {}".format(inputfile, geometrytype, desc.shapeType))
        return inputfile


def OverlapCalc(HU_Input, Bloecke):
    """
    :param HU_Input:
    :param Bloecke:
    :param OverlapThreshold:
    :return: City blocks and related buildings below given OverlapThreshold

    - calculates Overlap ratio of sum of building footprints area to area of city block for each block
    - returns buildings and city blocks below OverlapThreshold
    """

    # locale variables
    HU_input_SP = mem("HU_input_SP")
    HU_Rectangle = mem("HU_Rectangle")
    HU_Rect_Buff = mem("HU_Rect_Buff")
    Bloeck_Merge = tmp("Bloeck_Merge.shp")
    Bloeck_Merge_Join = tmp("Bloeck_Merge_Join.shp")
    Bloeck_Merge_Join_Dissolve = tmp("Bloeck_Merge_Join_Dissolve.shp")

    # order buildings to the blocks
    MergeList = arcpy.ListFeatureClasses('Block_*')
    DelName(MergeList)
    if len(arcpy.ListFields(Bloecke, 'NAME')) > 0:
        arcpy.DeleteField_management(Bloecke, 'NAME')
    arcpy.AddField_management(Bloecke, 'NAME', 'TEXT')
    arcpy.CalculateField_management(Bloecke, "NAME", '"Block_" + str(!FID!)', "PYTHON_9.3", "#")
    arcpy.management.MultipartToSinglepart(HU_Input, HU_input_SP)
    Bloecke_FL = arcpy.MakeFeatureLayer_management(Bloecke)
    arcpy.analysis.Split(HU_input_SP, Bloecke_FL, "NAME", Workspace, None)

    MergeList = arcpy.ListFeatureClasses('Block_*')

    for i in MergeList:
        file = CheckFileType(i, "ShapeFile", "Polygon")
        if file != None:
            print " {} is None-Type and deleted.".format(file)
            arcpy.Delete_management(file)

    arcpy.management.Merge(MergeList, Bloeck_Merge)
    DelName(MergeList)
    arcpy.analysis.SpatialJoin(Bloeck_Merge, Bloecke, Bloeck_Merge_Join, "JOIN_ONE_TO_ONE", "KEEP_ALL", None,
                               "WITHIN", None, None)
    arcpy.management.Dissolve(Bloeck_Merge_Join, Bloeck_Merge_Join_Dissolve, "NAME", None, "MULTI_PART",
                              "DISSOLVE_LINES")

    # calculate area values of building footprints and blocks
    arcpy.AddField_management(Bloeck_Merge_Join_Dissolve, "AREA_BLK_M", "DOUBLE")
    arcpy.AddField_management(Bloecke, "SHAPE_AREA", "DOUBLE")
    arcpy.management.CalculateField(Bloeck_Merge_Join_Dissolve, "AREA_BLK_M", "!shape.geodesicArea@SQUAREMETERS!",
                                    "PYTHON_9.3", None)
    arcpy.management.CalculateField(Bloecke, "SHAPE_AREA", "!shape.geodesicArea@SQUAREMETERS!", "PYTHON_9.3", None)
    arcpy.management.JoinField(Bloecke, "NAME", Bloeck_Merge_Join_Dissolve, "NAME", "AREA_BLK_M")
    if len(arcpy.ListFields(Bloecke, 'OVERLAP')) == 0:
        arcpy.AddField_management(Bloecke, "OVERLAP", "DOUBLE", 5, "", "", "OVERLAP", "NULLABLE")

    # calculate overlap
    i = 0
    with arcpy.da.UpdateCursor(Bloecke, ["OVERLAP", "AREA_BLK_M", "SHAPE_AREA"]) as cursor18:
        for row in cursor18:
            i = i + 1
            try:
                row[0] = row[1] / row[2] * 100
                cursor18.updateRow(row)
            except:
                print ("Overlap skipped at row {}!!!!".format(i))
    del cursor18

    DelName([HU_input_SP, HU_Rectangle, HU_Rect_Buff, Bloeck_Merge, Bloeck_Merge_Join, "Blocks_FL", \
             Bloeck_Merge_Join_Dissolve])

    return Bloecke


def stat(inputdata):
    Shp_Area(inputdata)
    arcpy.analysis.Statistics(inputdata, "inputstattabl.dbf", "Shape_Area SUM")

    frequency, area = 0, 0

    with arcpy.da.SearchCursor("inputstattabl.dbf", ["FREQUENCY", "SUM_Shape_"]) as cursor1:
        for row in cursor1:
            frequency, area = row
    del cursor1

    area = round(area / 10000)
    DelName(["inputstattabl.dbf"])
    return frequency, area


def cls():
    os.system('cls' if os.name == 'nt' else 'clear')


def preparation(UGB, IBS, InputHU):
    SymDiff = tmp('SymDiff.shp')
    SymDiffSP = tmp('SymDiffSP.shp')

    UGB_FL = arcpy.MakeFeatureLayer_management(UGB)
    arcpy.analysis.SymDiff(UGB, IBS, SymDiff, "ALL", None)

    if len(arcpy.ListFields(SymDiff, "INNEN")) > 0:
        arcpy.DeleteField_management(SymDiff, "INNEN")
    arcpy.AddField_management(SymDiff, "INNEN", "TEXT")
    arcpy.management.MultipartToSinglepart(SymDiff, SymDiffSP)
    SymDiffSP_FL = arcpy.MakeFeatureLayer_management(SymDiffSP)

    # Calc POS / NEG
    arcpy.management.SelectLayerByLocation(SymDiffSP_FL, "WITHIN", UGB_FL, None, "NEW_SELECTION", "NOT_INVERT")
    arcpy.management.CalculateField(SymDiffSP_FL, "INNEN", "'POS'", "PYTHON_9.3")
    arcpy.SelectLayerByAttribute_management(SymDiffSP_FL, "SWITCH_SELECTION")
    arcpy.management.CalculateField(SymDiffSP_FL, "INNEN", "'NEG'", "PYTHON_9.3")

    Shp_Length(SymDiffSP)
    Shp_Area(SymDiffSP)


    # NAME
    if len(arcpy.ListFields(SymDiffSP, "NAME")) > 0:
        arcpy.DeleteField_management(SymDiffSP, "NAME")
    arcpy.AddField_management(SymDiffSP, "NAME", "TEXT")
    arcpy.CalculateField_management(SymDiffSP, 'NAME', '"Block_" + str(!FID!)', 'PYTHON_9.3', '#')

    # OVERLAP

    SymDiffOC = OverlapCalc(InputHU, SymDiffSP)

    arcpy.RepairGeometry_management(SymDiffOC)


    return SymDiffOC


def Divide_POS_NEG(SymDiff):
    SymDiff_POS = tmp('SymDiff_POS.shp')
    SymDiff_NEG = tmp('SymDiff_NEG.shp')

    SymDiff_FL = arcpy.MakeFeatureLayer_management(SymDiff)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiff_FL, "NEW_SELECTION", "INNEN = 'POS'")
    arcpy.CopyFeatures_management(Sel, SymDiff_POS)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiff_FL, "NEW_SELECTION", "INNEN = 'NEG'")
    arcpy.CopyFeatures_management(Sel, SymDiff_NEG)

    return SymDiff_POS, SymDiff_NEG


def Class_IndCom(SymDiff, Nutzungen, GOT, INNEN):
    '''

    :param SymDiff:
    :param Nutzungen:
    :param GOT:
    :return:
    '''

    IndGew = mem('IndGew')
    SymDiffIndGew = mem('SymDiffIndGew')
    SymDiff_P1 = tmp("A_{}_IndCom.shp".format(INNEN))
    SymDiff_P1_IS = mem('SymDiff_{}_IS'.format(INNEN))
    SymDiff_P1_IS_DISS = mem('SymDiff_{}_IS_DISS'.format(INNEN))
    SymDiff_IS_Sel = mem('SymDiff_IS_Sel_{}'.format(INNEN))
    SymDiff_join = mem('SymDiff_join_{}'.format(INNEN))
    SymDiff_P1_join = mem('SymDiff_P1_join_{}'.format(INNEN))

    Nutzungen_FL = arcpy.MakeFeatureLayer_management(Nutzungen)
    Sel = arcpy.management.SelectLayerByAttribute(Nutzungen_FL, "NEW_SELECTION",
                                                  "OBJART_TXT = 'AX_IndustrieUndGewerbeflaeche'")
    arcpy.CopyFeatures_management(Sel, IndGew)
    IndGew_FL = arcpy.MakeFeatureLayer_management(IndGew)
    SymDiff_FL = arcpy.MakeFeatureLayer_management(SymDiff)
    Sel = arcpy.management.SelectLayerByLocation(SymDiff_FL, "INTERSECT", IndGew_FL, None, "NEW_SELECTION",
                                                 "NOT_INVERT")
    arcpy.CopyFeatures_management(Sel, SymDiffIndGew)
    Shp_Area(SymDiffIndGew)
    SymDiffIndGew_FL = arcpy.MakeFeatureLayer_management(SymDiffIndGew)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiffIndGew_FL, "NEW_SELECTION",
                                                  "INNEN = '{}' And  Shape_Area > 10000 And OVERLAP > {}".format(INNEN,
                                                                                                                 GOT))
    arcpy.CopyFeatures_management(Sel, SymDiff_P1)
    Shp_Area(SymDiff_P1)

    arcpy.analysis.Intersect([SymDiff_P1, IndGew], SymDiff_P1_IS, "ALL", None, "INPUT")
    if len(arcpy.ListFields(SymDiff_P1, "ORIG_FID")) > 0:
        arcpy.DeleteField_management(SymDiff_P1, "ORIG_FID")
    if len(arcpy.ListFields(SymDiff_P1_IS, "ORIG_FID")) > 0:
        arcpy.DeleteField_management(SymDiff_P1_IS, "ORIG_FID")
    arcpy.analysis.SpatialJoin(SymDiff_P1_IS, SymDiff_P1, SymDiff_P1_join)
    arcpy.management.Dissolve(SymDiff_P1_join, SymDiff_P1_IS_DISS, "FID_A_{}_IndCom".format(INNEN), None,
                              "MULTI_PART", "DISSOLVE_LINES")

    Shp_Area(SymDiff_P1_IS_DISS)
    arcpy.AlterField_management(SymDiff_P1_IS_DISS, "Shape_Area", "Ind_Area")
    arcpy.analysis.SpatialJoin(SymDiff_P1, SymDiff_P1_IS_DISS, SymDiff_join)
    if len(arcpy.ListFields(SymDiff_join, "RATIO")) == 0:
        arcpy.management.AddField(SymDiff_join, "RATIO", "DOUBLE")
    arcpy.management.CalculateField(SymDiff_join, "RATIO", "!Ind_Area! / !Shape_Area! * 100", "PYTHON_9.3", None)
    SymDiff_join_FL = arcpy.MakeFeatureLayer_management(SymDiff_join)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiff_join_FL, "NEW_SELECTION", "RATIO < 50")
    arcpy.CopyFeatures_management(Sel, SymDiff_IS_Sel)
    SymDiff_IS_Sel_FL = arcpy.MakeFeatureLayer_management(SymDiff_IS_Sel)
    SymDiff_P1_FL = arcpy.MakeFeatureLayer_management(SymDiff_P1)

    arcpy.management.SelectLayerByLocation(SymDiff_P1_FL, "ARE_IDENTICAL_TO", SymDiff_IS_Sel_FL, None, "NEW_SELECTION",
                                           "NOT_INVERT")
    arcpy.CopyFeatures_management(SymDiff_P1, "SymDiff_P1_{}_0.shp".format(INNEN)) #############
    arcpy.CopyFeatures_management(SymDiff_IS_Sel_FL, "SymDiff_IS_Sel_FL_{}_0.shp".format(INNEN))  #############
    arcpy.DeleteFeatures_management(SymDiff_P1_FL)
    arcpy.CopyFeatures_management(SymDiff_FL, "SymDiff_FL_{}_0.shp".format(INNEN))  #############
    arcpy.management.SelectLayerByLocation(SymDiff_FL, "ARE_IDENTICAL_TO", SymDiff_P1_FL, None, "NEW_SELECTION",
                                           "NOT_INVERT")
    arcpy.DeleteFeatures_management(SymDiff_FL)
    arcpy.CopyFeatures_management(SymDiff_FL, "SymDiff_FL_{}_1.shp".format(INNEN))  #############

    return SymDiff_FL, SymDiff_P1


def Class_Holes(SymDiff, IBS, INNEN):

    SymDiff_N3 = tmp("A_{}_Holes.shp".format(INNEN))
    UGBDiss = mem('UGBDiss')
    UGB_holes = mem('UGBHoles')
    IBS_Diss = mem('IBS_Diss')
    IBS_holes = mem('IBSHoles')

    # identify all holes in UGB
    arcpy.management.Dissolve(UGB, UGBDiss, None, None, "SINGLE_PART")
    UGB_FL = arcpy.MakeFeatureLayer_management(UGBDiss)
    arcpy.management.FeatureToLine(UGBDiss, mem("UGBLine"), None, "ATTRIBUTES")
    arcpy.management.FeatureToPolygon(mem("UGBLine"), mem("UGBLinePoly"), None, "ATTRIBUTES", None)
    UGBLinePoly_FL = arcpy.MakeFeatureLayer_management(mem("UGBLinePoly"))
    Sel = arcpy.management.SelectLayerByLocation(UGBLinePoly_FL, "ARE_IDENTICAL_TO", UGB_FL, None, "NEW_SELECTION",
                                                 "INVERT")
    arcpy.CopyFeatures_management(Sel, UGB_holes)
    # identify all holes in IBS
    arcpy.management.Dissolve(IBS, IBS_Diss, None, None, "SINGLE_PART")
    IBS_FL = arcpy.MakeFeatureLayer_management(IBS_Diss)
    arcpy.management.FeatureToLine(IBS_Diss, mem("IBSLine"), None, "ATTRIBUTES")
    arcpy.management.FeatureToPolygon(mem("IBSLine"), mem("IBSLinePoly"), None, "ATTRIBUTES", None)
    IBSLinePoly_FL = arcpy.MakeFeatureLayer_management(mem("IBSLinePoly"))
    Sel = arcpy.management.SelectLayerByLocation(IBSLinePoly_FL, "ARE_IDENTICAL_TO", IBS_FL, None, "NEW_SELECTION",
                                                 "INVERT")
    arcpy.CopyFeatures_management(Sel, IBS_holes)

    SymDiff_FL = arcpy.MakeFeatureLayer_management(SymDiff)
    # case select
    if INNEN == 'POS':
        IBS_holes_FL = arcpy.MakeFeatureLayer_management(IBS_holes)
        Sel = arcpy.management.SelectLayerByLocation(SymDiff_FL, "ARE_IDENTICAL_TO", IBS_holes_FL, None,
                                                     "NEW_SELECTION", "NOT_INVERT")
        arcpy.CopyFeatures_management(Sel, mem("SymDiffHoles"))
        SymDiffHoles_FL = arcpy.MakeFeatureLayer_management(mem("SymDiffHoles"))
        Sel = arcpy.management.SelectLayerByAttribute(SymDiffHoles_FL, "NEW_SELECTION",
                                                      "OVERLAP < {} AND INNEN = '{}'".format(LBC, INNEN))

    if INNEN == 'NEG':
        UGBHoles_FL = arcpy.MakeFeatureLayer_management(UGB_holes)
        Sel = arcpy.management.SelectLayerByLocation(SymDiff_FL, "ARE_IDENTICAL_TO", UGBHoles_FL, None,
                                                     "NEW_SELECTION", "NOT_INVERT")
        arcpy.CopyFeatures_management(Sel, mem("SymDiffHoles"))
        SymDiffHoles_FL = arcpy.MakeFeatureLayer_management(mem("SymDiffHoles"))
        Sel = arcpy.management.SelectLayerByAttribute(SymDiffHoles_FL, "NEW_SELECTION",
                                                      "OVERLAP < {} AND INNEN = '{}'".format(LBC, INNEN))

    # create output files
    arcpy.CopyFeatures_management(Sel, SymDiff_N3)

    Shp_Area(SymDiff_N3)
    SymDiff_N3_FL = arcpy.MakeFeatureLayer_management(SymDiff_N3)
    arcpy.management.SelectLayerByLocation(SymDiff_FL, "ARE_IDENTICAL_TO", SymDiff_N3_FL, None,
                                           "NEW_SELECTION",
                                           "NOT_INVERT")
    arcpy.DeleteFeatures_management(SymDiff_FL)


    return SymDiff_FL, SymDiff_N3


def Class_SettBody(SymDiff, IBS, UGB, INNEN):
    """
    IBS where is no UGB
    """

    SymDiff_N7 = tmp("A_{}_SettBody.shp".format(INNEN))

    IBS_FL = arcpy.MakeFeatureLayer_management(IBS)
    UGB_FL = arcpy.MakeFeatureLayer_management(UGB)
    SymDiff_FL = arcpy.MakeFeatureLayer_management(SymDiff)

    if INNEN == 'POS':
        Sel = arcpy.management.SelectLayerByLocation(UGB_FL, "INTERSECT", IBS_FL, None, "NEW_SELECTION", "INVERT")
    else:
        Sel = arcpy.management.SelectLayerByLocation(IBS_FL, "INTERSECT", UGB_FL, None, "NEW_SELECTION", "INVERT")

    arcpy.CopyFeatures_management(Sel, SymDiff_N7)
    SymDiff_N7_FL = arcpy.MakeFeatureLayer_management(SymDiff_N7)
    arcpy.management.SelectLayerByLocation(SymDiff_N7_FL, "INTERSECT", SymDiff_FL , None,
                                           "NEW_SELECTION",
                                           "INVERT")
    arcpy.DeleteFeatures_management(SymDiff_N7_FL)
    Shp_Area(SymDiff_N7)
    SymDiff_N7_FL = arcpy.MakeFeatureLayer_management(SymDiff_N7)
    arcpy.management.SelectLayerByLocation(SymDiff_FL, "ARE_IDENTICAL_TO", SymDiff_N7_FL, None,
                                           "NEW_SELECTION",
                                           "NOT_INVERT")
    arcpy.DeleteFeatures_management(SymDiff_FL)

    return SymDiff_FL, SymDiff_N7


def Class_Resid(SymDiff, Nutzungen, GOT, INNEN):

    SymDiffPos = mem('SymDiffPos')
    WohnGemFunk = tmp("WohnGemFunk.shp")
    SymDiffPos_Wohn = tmp(('SymDiffPos_Wohn.shp'))

    SymDiff_P2 = tmp("A_{}_Resid.shp".format(INNEN))
    SymDiff_P2_IS = mem('SymDiff_{}1_IS'.format(INNEN))
    SymDiff_P2_IS_DISS = mem('SymDiff_{}1_IS_DISS'.format(INNEN))
    SymDiff_IS_Sel = mem('SymDiff_IS_Sel_{}'.format(INNEN))
    SymDiff_join = mem('SymDiff_join_{}'.format(INNEN))
    SymDiff_P2_join = mem('SymDiff_P2_join_{}'.format(INNEN))

    SymDiff_FL = arcpy.MakeFeatureLayer_management(SymDiff)
    Nutzungen_FL = arcpy.MakeFeatureLayer_management(Nutzungen)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiff_FL, "NEW_SELECTION", "INNEN = '{}'".format(INNEN))
    arcpy.CopyFeatures_management(Sel, SymDiffPos)
    Shp_Area(SymDiffPos)
    SymDiffPos_FL = arcpy.MakeFeatureLayer_management(SymDiffPos)
    Sel = arcpy.management.SelectLayerByAttribute(Nutzungen_FL, "NEW_SELECTION",
                                                  "OBJART_TXT = 'AX_Wohnbauflaeche' Or OBJART_TXT = 'AX_FlaecheGemischterNutzung' Or OBJART_TXT = 'AX_FlaecheBesondererFunktionalerPraegung'")
    arcpy.CopyFeatures_management(Sel, WohnGemFunk)
    WohnGemFunk_FL = arcpy.MakeFeatureLayer_management(WohnGemFunk)
    Sel = arcpy.management.SelectLayerByLocation(SymDiffPos_FL, "INTERSECT", WohnGemFunk_FL, None, "NEW_SELECTION",
                                                 "NOT_INVERT")
    arcpy.CopyFeatures_management(Sel, SymDiffPos_Wohn)
    Shp_Area(SymDiffPos_Wohn)

    SymDiffPos_Wohn_FL = arcpy.MakeFeatureLayer_management(SymDiffPos_Wohn)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiffPos_Wohn_FL, "NEW_SELECTION",
                                                  "OVERLAP > {} And Shape_Area > 10000".format(GOT))
    arcpy.CopyFeatures_management(Sel, SymDiff_P2)

    arcpy.analysis.Intersect([SymDiff_P2, WohnGemFunk], SymDiff_P2_IS, "ALL", None, "INPUT")

    arcpy.analysis.SpatialJoin(SymDiff_P2_IS, SymDiff_P2, SymDiff_P2_join)
    arcpy.management.Dissolve(SymDiff_P2_join, SymDiff_P2_IS_DISS, "ORIG_FID", None, "MULTI_PART", "DISSOLVE_LINES")

    Shp_Area(SymDiff_P2_IS_DISS)
    arcpy.AlterField_management(SymDiff_P2_IS_DISS, "Shape_Area", "Res_Area")
    arcpy.analysis.SpatialJoin(SymDiff_P2, SymDiff_P2_IS_DISS, SymDiff_join,
                               "JOIN_ONE_TO_ONE", "KEEP_ALL", None, "INTERSECT", None, '')
    if len(arcpy.ListFields(SymDiff_join, "RATIO")) == 0:
        arcpy.management.AddField(SymDiff_join, "RATIO", "DOUBLE")
    arcpy.management.CalculateField(SymDiff_join, "RATIO", "!Res_Area! / !Shape_Area! * 100", "PYTHON_9.3", None)
    SymDiff_join_FL = arcpy.MakeFeatureLayer_management(SymDiff_join)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiff_join_FL, "NEW_SELECTION", "RATIO < 50")
    arcpy.CopyFeatures_management(Sel, SymDiff_IS_Sel)
    arcpy.CopyFeatures_management(Sel, SymDiff_IS_Sel)
    SymDiff_IS_Sel_FL = arcpy.MakeFeatureLayer_management(SymDiff_IS_Sel)
    SymDiff_P2_FL = arcpy.MakeFeatureLayer_management(SymDiff_P2)
    arcpy.management.SelectLayerByLocation(SymDiff_P2_FL, "ARE_IDENTICAL_TO", SymDiff_IS_Sel_FL, None, "NEW_SELECTION",
                                           "NOT_INVERT")
    arcpy.DeleteFeatures_management(SymDiff_P2_FL)
    arcpy.management.SelectLayerByLocation(SymDiff_FL, "ARE_IDENTICAL_TO", SymDiff_P2_FL, None,
                                           "NEW_SELECTION", "NOT_INVERT")
    arcpy.DeleteFeatures_management(SymDiff_FL)

    DelName([WohnGemFunk, SymDiffPos_Wohn])

    return SymDiff_FL, SymDiff_P2


def Class_BdgEdge(SymDiff, INNEN):
    """Build up areas O > 15 % """

    SymDiff_P4 = tmp("A_{}_BdgEdge.shp".format(INNEN))

    SymDiff_FL = arcpy.MakeFeatureLayer_management(SymDiff)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiff_FL, "NEW_SELECTION",
                                                  "OVERLAP > {} And INNEN = '{}'".format(GOT, INNEN))
    arcpy.CopyFeatures_management(Sel, SymDiff_P4)
    arcpy.DeleteFeatures_management(SymDiff_FL)

    return SymDiff_FL, SymDiff_P4


def Class_EmptyAreas(SymDiff, INNEN):
    """emty positve areas, O < 3% """

    SymDiff_P5 = tmp("A_{}_EmptyArea.shp".format(INNEN))

    SymDiff_FL = arcpy.MakeFeatureLayer_management(SymDiff)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiff_FL, "NEW_SELECTION", "OVERLAP <= {} And INNEN = '{}'".format(LBC, INNEN))
    arcpy.CopyFeatures_management(Sel, SymDiff_P5)
    arcpy.DeleteFeatures_management(SymDiff_FL)

    return SymDiff_FL, SymDiff_P5


def Class_LowDensBdgGrp(SymDiff, GOT, INNEN):

    SymDiff_P6 = tmp("A_{}_LowDensBdgGrp.shp".format(INNEN))

    SymDiff_FL = arcpy.MakeFeatureLayer_management(SymDiff)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiff_FL, "NEW_SELECTION",
                                                  "OVERLAP > {} And OVERLAP <= {} And INNEN = '{}'".format(LBC, GOT, INNEN))
    arcpy.CopyFeatures_management(Sel, SymDiff_P6)
    arcpy.DeleteFeatures_management(SymDiff_FL)

    return SymDiff_FL, SymDiff_P6


def Class_LargeEmptyAreas(SymDiff, INNEN):


    SymDiff_N5 = tmp("A_{}_LargeEmptyAreas.shp".format(INNEN))

    SymDiff_FL = arcpy.MakeFeatureLayer_management(SymDiff)
    Sel = arcpy.management.SelectLayerByAttribute(SymDiff_FL, "NEW_SELECTION",
                                                  "OVERLAP <= {} And Shape_Area > 10000 And INNEN = '{}'".format(LBC, INNEN))
    arcpy.CopyFeatures_management(Sel, SymDiff_N5)
    arcpy.DeleteFeatures_management(SymDiff_FL)

    return SymDiff_FL, SymDiff_N5


def main():

    IntersFreature = mem('IntersFreature')
    IBSSP = mem('IBSSP')
    SymDiff_Rest = tmp('SymDiff_Rest.shp')
    SymDiff_Copy = tmp('SymDiff_Copy.shp')

    DelName(["Tmp"])
    arcpy.CreateFolder_management(Workspace, "Tmp")

    Shp_Area(UGB)
    Shp_Area(IBS)

    SymDiff = preparation(UGB, IBS, InputHU)

    arcpy.CopyFeatures_management(SymDiff, SymDiff_Copy)

    arcpy.Intersect_analysis([UGB, IBS], IntersFreature)
    arcpy.management.MultipartToSinglepart(IBS, IBSSP)

    Shp_Area(IntersFreature)
    anz_IntersF, area_IntersF = stat(IntersFreature)
    anz_UGB, area_UGB = stat(UGB)
    anz_IBS, area_IBS = stat(IBSSP)

    SymDiff_POS, SymDiff_NEG = Divide_POS_NEG(SymDiff_Copy)
    anz_POS, area_POS = stat(SymDiff_POS)
    anz_NEG, area_NEG = stat(SymDiff_NEG)

    sum_area = area_POS + area_NEG + area_IntersF

    share_IntersF = round(area_IntersF / sum_area * 100, 1)

    row_format = "{:<20} {:>6} {:>6} {:>6}"

    data = [
        ["Toatl area:", "", sum_area, ""],
        ["IBS Freq. Area:", anz_IBS, area_IBS, ""],
        ["UGB Freq. Area:", anz_UGB, area_UGB, ""],
        ["Slice Freq. Area:", anz_IntersF, area_IntersF, share_IntersF],
        ["", "", "", ""]]

    print (row_format.format("", "Frequ.", "Area", "Share"))
    for row in data:
        print(row_format.format(*row))

    # delete patches < 250 m2
    SymDiff_FL = arcpy.MakeFeatureLayer_management(SymDiff)
    arcpy.management.SelectLayerByAttribute(SymDiff_FL, "NEW_SELECTION", "Shape_Area < 250")
    arcpy.DeleteFeatures_management(SymDiff_FL)

    SymDiff, SymDiff_P1 = Class_IndCom(SymDiff, Nutzungen, GOT, 'POS')

    SymDiff, SymDiff_N_IndCom = Class_IndCom(SymDiff, Nutzungen, GOT, 'NEG')

    SymDiff, SymDiff_N1 = Class_Resid(SymDiff, Nutzungen, GOT, 'NEG')

    SymDiff, SymDiff_N3 = Class_Holes(SymDiff, IBS, 'NEG')
    SymDiff, SymDiff_N7 = Class_SettBody(SymDiff, IBS, UGB, 'NEG')

    SymDiff, SymDiff_P2 = Class_Resid(SymDiff, Nutzungen, GOT, 'POS')

    SymDiff, SymDiff_P_holes = Class_Holes(SymDiff, IBS, 'POS')

    SymDiff, SymDiff_P_SettBody = Class_SettBody(SymDiff, IBS, UGB, 'POS')

    SymDiff, SymDiff_P_LargeEmpty = Class_LargeEmptyAreas(SymDiff, 'POS')

    SymDiff, SymDiff_P4 = Class_BdgEdge(SymDiff, 'POS')

    SymDiff, SymDiff_P5 = Class_EmptyAreas(SymDiff, 'POS')

    SymDiff, SymDiff_P6 = Class_LowDensBdgGrp(SymDiff, GOT, 'POS')

    SymDiff, SymDiff_N5 = Class_LargeEmptyAreas(SymDiff, 'NEG')

    SymDiff, SymDiff_N8 = Class_BdgEdge(SymDiff, 'NEG')

    SymDiff, SymDiff_N4 = Class_EmptyAreas(SymDiff, 'NEG')

    SymDiff, SymDiff_N2 = Class_LowDensBdgGrp(SymDiff, GOT, 'NEG')


    arcpy.CopyFeatures_management(SymDiff, SymDiff_Rest)
    result = arcpy.GetCount_management(SymDiff_Rest)
    Anz = int(result.getOutput(0))

    SymDiff_POS, SymDiff_NEG = Divide_POS_NEG(SymDiff_Copy)
    anz_POS, area_POS = stat(SymDiff_POS)
    anz_P1, area_P1 = stat(SymDiff_P1)
    anz_P2, area_P2 = stat(SymDiff_P2)
    anz_P4, area_P4 = stat(SymDiff_P4)
    anz_P5, area_P5 = stat(SymDiff_P5)
    anz_P6, area_P6 = stat(SymDiff_P6)
    anz_P_holes, area_P_holes = stat(SymDiff_P_holes)

    anz_P_LargeEmpty, area_P_LargeEmpty = stat(SymDiff_P_LargeEmpty)
    anz_P_SettBody, area_P_SettBody = stat(SymDiff_P_SettBody)

    sum_anz_P = anz_P1 + anz_P2 + anz_P4 + anz_P5 + anz_P6 + anz_P_LargeEmpty + anz_P_SettBody + anz_P_holes
    sum_area_P = area_P1 + area_P2 + area_P4 + area_P5 + area_P6 + area_P_LargeEmpty + area_P_SettBody + area_P_holes

    anz_NEG, area_NEG = stat(SymDiff_NEG)
    anz_N_IndCom, area_N_IndCom = stat(SymDiff_N_IndCom)
    anz_N1, area_N1 = stat(SymDiff_N1)
    anz_N2, area_N2 = stat(SymDiff_N2)
    anz_N3, area_N3 = stat(SymDiff_N3)
    anz_N4, area_N4 = stat(SymDiff_N4)
    anz_N5, area_N5 = stat(SymDiff_N5)
    anz_N7, area_N7 = stat(SymDiff_N7)
    anz_N8, area_N8 = stat(SymDiff_N8)

    sum_anz_N = anz_N1 + anz_N2 + anz_N3 + anz_N4 + anz_N5 + anz_N7 + anz_N8 + anz_N_IndCom
    sum_area_N = area_N1 + area_N2 + area_N3 + area_N4 + area_N5 + area_N7 + area_N8 + area_N_IndCom

    share_P1 = round(area_P1 / sum_area * 100, 1)
    share_P2 = round(area_P2 / sum_area * 100, 1)
    share_P5 = round(area_P5 / sum_area * 100, 1)
    share_P4 = round(area_P4 / sum_area * 100, 1)
    share_P6 = round(area_P6 / sum_area * 100, 1)
    share_P_LargeEmpty = round(area_P_LargeEmpty / sum_area * 100, 1)
    share_P_SettBody = round(area_P_SettBody / sum_area * 100, 1)
    share_P_holes = round(area_P_holes / sum_area * 100, 1)

    share_N_IndCom = round(area_N_IndCom / sum_area * 100, 1)
    share_N1 = round(area_N1 / sum_area * 100, 1)
    share_N2 = round(area_N2 / sum_area * 100, 1)
    share_N3 = round(area_N3 / sum_area * 100, 1)
    share_N4 = round(area_N4 / sum_area * 100, 1)
    share_N5 = round(area_N5 / sum_area * 100, 1)
    share_N7 = round(area_N7 / sum_area * 100, 1)
    share_N8 = round(area_N8 / sum_area * 100, 1)

    sum_share_P = share_P1 + share_P2 + share_P5 + share_P4 + share_P6 + share_P_LargeEmpty + share_P_SettBody + share_P_holes
    sum_share_N = share_N1 + share_N2 + share_N3 + share_N4 + share_N5 + share_N7 + share_N8 + share_N_IndCom

    data = [
        ["Class P_IndCom:", anz_P1, area_P1, share_P1],
        ["Class P_Resid:", anz_P2, area_P2, share_P2],
        ["Class P_BdgEdg:", anz_P4, area_P4, share_P4],
        ["Class P_LowDens:", anz_P6, area_P6, share_P6],
        ["Class P_LargeEmpty:", anz_P_LargeEmpty, area_P_LargeEmpty, share_P_LargeEmpty],
        ["Class P_EmptyArea:", anz_P5, area_P5, share_P5],
        ["Class P_SettBody:", anz_P_SettBody, area_P_SettBody, share_P_SettBody],
        ["Class P_Holes:", anz_P_holes, area_P_holes, share_P_holes],
        ["Sum:", sum_anz_P, sum_area_P, sum_share_P],
        ["", "", "", ""],
        ["Class N_IndCom:", anz_N_IndCom, area_N_IndCom, share_N_IndCom],
        ["Class N_Resid:", anz_N1, area_N1, share_N1],
        ["Class N_BdgGrp:", anz_N8, area_N8, share_N8],
        ["Class N_LowDens:", anz_N2, area_N2, share_N2],
        ["Class N_LargeEmpty:", anz_N5, area_N5, share_N5],
        ["Class N_EmptyAreas:", anz_N4, area_N4, share_N4],
        ["Class N_SettBoddy:", anz_N7, area_N7, share_N7],
        ["Class N_Holes:", anz_N3, area_N3, share_N3],
        ["Sum:", sum_anz_N, sum_area_N, sum_share_N],
        ["", "", "", ""],
        ["Unclassifid:", anz_POS + anz_NEG - sum_anz_P - sum_anz_N, sum_area - sum_area_P - sum_area_N - area_IntersF,
         100 - sum_share_N - sum_share_P - share_IntersF]]

    print (row_format.format("Class", "Frequ.", "Area", "Share"))
    for row in data:
        print(row_format.format(*row))

    if Anz > 0:
        print ("There are " + str(Anz) + "features left!")

    DelName([SymDiff])

if __name__ == '__main__':
    main()
