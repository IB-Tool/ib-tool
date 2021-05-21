# IB-Tool

Toolset for the delineation of settlements on the basis building footprints, road network and land use data

This tool was developed as part of a Phd thesis. For further information, please see the following publication: 
Harig, O.; Hecht, R.; Burghardt, D.; Meinel, G. Automatic Delineation of Urban Growth Boundaries Based on Topographic Data Using Germany as a Case Study. ISPRS Int. J. Geo-Inf. 2021, 10, 353. https://doi.org/10.3390/ijgi10050353 


## 1. Prerequisites

- Install the necessary packages via pip or the package maganger of your IDE.
- You need a running ArcGIS version (ArcMap or ArcGIS Prso) on the system on which the code is to be executed. Furthermore, a [Spatial Analyst licence](https://desktop.arcgis.com/en/arcmap/latest/extensions/spatial-analyst/what-is-the-spatial-analyst-extension.htm) is required.


### Python and Python packages 

In addition, the following Python packages have to be installed:

- Python == 2.7
- traceback
- networkx
- numpy
- scipy
- statistics
- msvcrt
- re

The script also uses the following standard modules:

- os
- time
- datetime
- operator 

## 2. Installing

The tool consists of three independent scripts:

1. partitioning.pyt - To create partitioning from building footprints
2. IB-Tool2.py, IB-Tool2_Config.txt, IB-Tool2_Filter.txt - To create the delineation
3. Error_Classification.py - To assess the boundary on the basis of expert delineation

It is recommended to create a separate directory for each of the three parts.

## 3. Application



### a) Partitioning

To be able to create a delineation, a partitioning must be created. To do this, the file partitioning.pyt from Python-Toolbox must be opened in ArcGIS. The building data can then be selected via the dialogue.
These files must be located in a geodatabase. The working directory must also be a geodatabase. It is advised not to change the default parameters.

If other data is used for partitioning, care must be taken that it contains a "NAME" field (text). The field "NAME" must contain an individual value for each feature, e.g. PART_1, PART_2, ...

### b) Creating settlement boundary

The script is designed to use ATKIS data as input data. If other data are used, problems may occur due to different field names/data formatting within the data.

Copy building floor plans, road network data, auxiliary data, partitioning data and the three script files(IB-Tool2.py, IB-Tool2_Config.txt, IB-Tool2_Filter.txt) into a working directory.

### c) Data preparation

#### Building Footprint

The building ground plans must be available as a polygon shape file. The data must not contain a "NAME" field.
The file must be renamed to A_HU.shp.

#### Road network data

The data must be available as a polyline shape file.
The file must be renamed to A_RN.shp.

#### Auxiliary data

The data must be available as a polyline shape file.
The file must be renamed to A_AUX.shp.

#### Partitioning data

The data must be available as a polygon shape file. The data must contain a "NAME" field.
The file must be renamed to A_PART.shp.

### d) Configuration/Parameterisation

#### Filters

Building filters are set in the IB-Tool2_Filter.txt file.

The identifier (digit) is decisive for the filtering. It is composed of the digit for the identifier (Kennung) and the digit for the value (Wert) according to [ATKIS-OK](http://www.adv-online.de/icc/extdeu/nav/a63/binarywriterservlet%3FimgUid%3D9201016e-7efa-8461-e336-b6951fa2e0c9%26uBasVariant%3D11111111-1111-1111-1111-111111111111).
Entries after the comma are used for commenting.

```
#Filter positive (AX_Gebaedefunktionen: Kennung, Wert, Bezeichner) jeweils ohne Leerzeichen eingeben
31001_1000, Wohngeb
31001_1010, Wohnhaus
31001_1100, GemischtesWohnen
31001_1120, WohnenHandelDienst
31001_1121, WohnVerw
31001_1122, WohnBuero
31001_1123, WohnGesch
31001_2052, Einkaufszentrum
31001_2010, Handel und Dienstleistungen
31001_2050, Geschäftsgebäude
31001_2020, Bürogebäude
31001_3021, Schule
31001_3000, öffentliche Zwecke
31001_3017, Kreisverwaltung
31001_3018, Bezirksregierung

#Filter negative (AX_Gebaedefunktionen: Kennung, Wert, Bezeichner) jeweils ohne Leerzeichen eingeben
31001_1310, Freizeit
31001_2600, Entsorgung
31001_2720, GebLandForst
31001_2720, GebLandForst
31001_2721, Scheune
31001_2723, Schuppen
31001_2724, Stall
31001_2726, StallSchuppen
31001_2727, Stall
31001_2740, Treibhaus
31001_2741, Treibhaus
31001_2742, Treibhaus
31001_2140, Vorratshaltung
51003_1201, Silo
31001_2143, Lager
51002_1215, Biogas
31001_3200, Erholung
31001_2463, Garage
31001_2523, Umfomer
31001_1312, Wochenendhaus
#End
```

#### Parameter

Parameters are set in IB-Tool2_Config.txt file.

```
01 Minimum Overlap of Inner Area Blocks (%): 18
02 Global Building Coverage Threshold (%): 0
03 Minimum Building Footprint Area (sqm): 56.8
04 Minimum Number of Buildings to create a Inner Zone: 20
05 Minimum Size of Inner Zone Patch (sqm): 10000
06 Maximum Size of Hole (sqm): 10000
07 Maximum Size of Gaps (sqm): 4900
08 Partition Start (default value = -1): -1
09 Partition End (default value = -1): -1
10 Partition List (default value = #):#
11 Log Level (default value = Warning):Debug
12 Spatial Reference: 25832
13 Path Common Workspace: D:\path\to\working\directory
14 Delete PartLog at start:True
```

1. Building Coverage Threshold for selecting dense blocks (default value = 18)
2. Global Building Coverage Threshold can be added if this has already been calculated (default value = 0) 
3. Minimum Building Footprint Area (default value = 56.8 sqm )
4. Minimum Number of Buildings to create a inner Zone/ UGB (default value = 20)
5. Minimum Size of Inner Zone Patch (default value = 10000 sqm)
6. Maximum Size of Hole (default value = 10000 sqm)
7. Maximum Size of Gaps (default value = 4900 sqm)
The values 8. 9, 10 are used for partioning purposes of the data set
8. Start value (OID or FID of feature in A_PART.shp) (default value = -1 = no filter)
9. End value (OID or FID of feature in A_PART.shp) (default value = -1 = no filter)
10. List of values according to "NAME" field of A_PART.shp (default value = #) for debugging or tests, no blanks
11. Set Log Level Alert, Warning or Debug (default value = Warning)
12. Set Spatial Reference by using a coordinate system's factory code (or authority code)
13. Set path to common workspace
14. Set switch for PartLog to True or False (default value = True)


### e) Comparison of expert delineations and calculated urban growth boundaries

Set imput data:

```
# ------ INPUT DATA ------

InputHU = Workspace + os.path.sep + "A_HU_P.shp" #  Building footprints (polygon, shape)
IBS = Workspace + os.path.sep + "A_IBS_P.shp" #  Expert delieation (polygon, shape)
UGB = Workspace + os.path.sep + "A_UGB_P5.shp"  #  Caculated boundary (polygon, shape)
Nutzungen = Workspace + os.path.sep + "A_Nutzungen_P.shp" #  land use geometry (polygon, shape)
sr = arcpy.SpatialReference(25833)
````

## 4. Authors

Oliver Harig: Conception, code implementation
Robert Hecht, Gotthardt Meinel, Dirk Burghardt: Conception

## 5. License

See [LICENSE.md](LICENSE.md) file for details.

