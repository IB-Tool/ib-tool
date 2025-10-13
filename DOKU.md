# Dokumentation des Tools zur automatisierten Siedlungsabgrenzung

*(nach Harig 2024, TU Dresden – Automatisierte Abgrenzung von Innenbereichen einschließlich Ergebnisevaluierung)*

---

## 1. Einleitung

Das vorliegende Tool dient der **automatisierten Abgrenzung von Siedlungs- bzw. Innenbereichen** auf Basis geometrischer und topographischer Objektdaten.
Es wurde im Rahmen der Dissertation von **Oliver Harig (2024)** entwickelt und zielt darauf ab, **einheitliche, datenbasierte und nachvollziehbare Siedlungsgrenzen** zu erzeugen, die den planerischen Innenbereich nach § 34 BauGB approximieren.

Die Methodik beruht auf der Kombination von:

* **Gebäudegrundrissen (HU-DE)** als Indikator realer Bebauung,
* **topographischen Linien (ATKIS-Straßennetz)** zur Partitionierung der Raumstruktur,
* **Hilfgesometrien (ATKIS-Wald, -Gewässer, -Moore u. Sumpf)**,
* **dichtebasierten Clustering-Algorithmen** zur Identifikation bebauter Bereiche,
* und **graphbasierten Aggregationsverfahren (Minimaler Spannbaum)** zur flächigen Zusammenführung.

Ziel ist eine **objektive und reproduzierbare Abgrenzung**, die zugleich planerisch interpretierbar bleibt.
Das Verfahren bildet eine Kernkomponente für **Siedlungsflächenmonitoring**, **Innenentwicklungspotenzialanalysen** und **Zersiedelungsstudien** (vgl. Harig 2024, Kap. 1.2, 1.4, 8.4).

---

## 2. Methodischer Ablauf des Verfahrens

Die Verarbeitung folgt in ihrer Grundstruktur dem in der Dissertation dargestellten Workflow (siehe Abb. 4.1–4.10 in Harig 2024).

1. **Partitionierung** des Untersuchungsgebiets in *Street Blocks* mittels Straßennetz.
2. **Erzeugung von City Blocks** durch Entfernen nicht geschlossener Partien (Abb. 4.6).
3. **Berechnung des Überbauungsgrads (UG)** pro Block auf Basis der Gebäudegrundflächen (Abb. 4.7).
4. **Filterung** kleiner oder irrelevanter Gebäude (z. B. Garagen, Nebengebäude).
5. **Identifikation dichter Blöcke** nach Schwellwerten (Kap. 4.5.3).
6. **Aggregation** mittels *Minimaler Spannbaum*-Algorithmus (Kap. 4.4.6).
7. **Nachbearbeitung**: Füllen von Baulücken und Verbinden benachbarter Siedlungen (Abb. 4.10).

Die im Tool konfigurierbaren **Parameter** greifen an verschiedenen Stellen dieses Prozesses ein und steuern dessen Empfindlichkeit, räumliche Kohärenz und Generalisierungsgrad.

---

## 3. Parameterbeschreibung

### 3.1 **`min_overlap_blocks` – Minimale Überlappung dichter Blöcke**

* **Typ:** Float
* **Funktion:** Steuert die Identifikation dicht bebauter Blöcke in `identify_dense_blocks()`.
* **Definition:** Gibt an, welcher Mindestanteil eines Blocks von Gebäudeflächen überlagert sein muss, um als „dicht bebaut“ klassifiziert zu werden.
* **Hintergrund:** Der Parameter operationalisiert den in der Dissertation beschriebenen **Überbauungsgrad (UG)**, der als Verhältnis von Gebäudegrundfläche zur Gesamtblockfläche definiert ist (vgl. Harig 2024, Kap. 4.4.3, S. 29 f.).
* **Wirkung:**

  * Bei **hohen Werten** (z. B. > 25) werden nur sehr kompakte, geschlossene Baublöcke als Siedlung erkannt (z.B. Stadtzentrum).
  * Bei **niedrigen Werten** (z. B. < 15) werden auch lockere Bebauungen (z. B. ländliche Bereiche) erfasst.
* **Empfohlener Bereich:** 18 – 22 (vgl. Kap. 4.5.3, empirische Kalibrierung).
* **Einfluss:** Dieser Schritt greift der Clusterung und Fangen der erzuegten Geometrien vor. Dadurch wird vorallem Rechenzeit gespart, da Baublöcke, die sicher als Innenbereich angesprochen werden können, in den nachfolgenden Schritten nicht betrachtet werden. Der Standardwert ist mit 18 relativ hoch gewählt, da so die Zahl der falsch-positive erkannten Baublöcke sehr niedrig ist. 

---

### 3.2 **`global_footprint_density` – Globale Gebäudedichte**

* **Typ:** Float
* **Funktion:** Wird als Referenzschwellwert für die Gesamtfläche verwendet, wenn keine lokalen Dichtewerte vorliegen (`calc_footprint_density()`).
* **Definition:** Mittelwert der Gebäudeflächendichte (Überbauungsgrad) aller Partitionen eines Untersuchungsgebiets.
* **Hintergrund:** Dient der Robustheit der Analyse in Gebieten mit geringer Bebauungsdichte oder unvollständigen Gebäudedaten (vgl. Harig 2024, Kap. 4.5.3).
* **Wirkung:** 

  * Niedrige Werte (< 10) → auch Gebäude in relativ lockere Bebauung werden aggregiert. 
  * Hohe Werte (> 20) → nur sehr verdichtete Siedlungskerne werden aggregiert.
* **Empfohlene Einstellung:** Automatische Berechnung (0 = auto); manuell bei stark heterogenen Untersuchungsgebieten setzen.
* **Bemerkung:** Der Wert wirkt sich per se nur auf die eher kleinen Siedlungen oder sehr dispersen Siedlungen aus. Er hat dann aber einen sehr großen Einfluss.

---

### 3.3 **`min_area` – Mindestgebäudefläche**

* **Typ:** Float
* **Funktion:** Filtert in `input_hu_filter()` zu kleine Gebäudeobjekte aus.
* **Definition:** Gebäude mit einer Grundfläche unterhalb dieses Werts werden für alle weiteren Berechnungen aus dem Datensatz gelöscht.
* **Ziel:** Ausschluss nicht siedlungsrelevanter Strukturen (z. B. Gartenhäuser, Carports, Nebengebäude).
* **Hintergrund:** Kleine Gebäude und Nebengebäude erfüllen nicht die Vorausetzungen, um einen Bebauungszusammenhang zu erzeugen. Sie werden daher vorab gelöscht.
* **Empfohlene Werte:** 50-60 m², je nach Örtlichkeit und regionalem Durchschnitt.
* **Wirkung:** Insbesondere am Siedlungsrand verläuft die Grenze dadurch direkt an den Wohngebäuden. Auf dem Grundstücken rückwärtige Gebäude werden so nicht berücksichtigt.

---

### 3.4 **`min_bdg_count` – Minimale Gebäudeanzahl**

* **Typ:** Integer
* **Funktion:** Definiert die Mindestanzahl von Gebäuden, ab der ein Cluster als gültiger Siedlungsbereich/Ortsteil gewertet wird.
* **Verwendung:** Wirkt sich erst in der Nachberarbeitung aus, indem Schwellwertabhängig Patsches gelöscht werden.
* **Hintergrund:** Laut empirischer Auswertung (Harig 2024, Kap. 4.5.6) und Expertenbefragung (Anhang A.1.3) gilt ein Ortsteil erst ab etwa **20 Gebäuden oder 1 ha** Fläche als eigenständiger Innenbereich.
* **Empfohlene Werte:** 15 – 25 Gebäude.
* **Wirkung:**

  * Verhindert Klassifikation isolierter Gebäudegruppen als Siedlung.
  * Stabilisiert Ergebnisse in ländlichen Räumen.

---

### 3.5 **`min_patch_size` – Minimale Patch-Größe**

* **Typ:** Float
* **Funktion:** Wird in `patch_remove()` zur Entfernung kleiner Siedlungsinseln eingesetzt.
* **Definition:** Unterhalb dieses Flächenwerts werden erzeugte Patches asl Siedlungsfläche verworfen.
* **Hintergrund:** Entspricht den empirisch ermittelten Mindestgrößen eines Ortsteils (Kap. 4.5.6).
* **Empfohlener Wert:** 10 000 m².
* **Wirkung:**

  * Kleine Fragmente werden entfernt → klarere Innenbereichsgrenzen.
  * Reduziert Übersegmentierung bei dünn besiedelten Regionen.
  * ABER: Flächenform wird nicht berücksichtigt!

---

### 3.6 **`max_hole_size` – Maximale Lochgröße**

* **Typ:** Float
* **Funktion:** In `gap_close()` zur Füllung unbebauter Flächen innerhalb zusammenhängender Siedlungsgeometrien.
* **Definition:** Gibt an, ab welcher Fläche Freiräume innerhalb eines Siedlungspolygons als *Außenbereich im Innenbereich* gelten bzw. gefüllt werden dürfen.
* **Hintergrund:** Entspricht planerischen Definitionen des Innenbereichs; 
* **Empfohlener Wert:** 10 000 m².
* **Wirkung:**

  * Kleinere Freiflächen (Gärten, Innenhöfe, Spielplätze, Wiesen, Parkplätze) werden integriert.
  * Verhindert übermäßige Fragmentierung der Siedlungsflächen.

---

### 3.7 **`max_gap_size` – Maximale Lückengröße**

* **Typ:** Float
* **Funktion:** In `gap_close()` schließt Lücken am Siedlungsrand.
* **Definition:** Gibt den Schwellwert an, bis zu welchem Lücken am Siedlungsrand geschlossen werden.
* **Hintergrund:** Baulücken zählen zum Innenbereich. Welche unbebaute Fläche noch als Baulücke gelten kann ist hoch umstritten und schwer zu fassen. Insbesondere sind die örtlichen Gegebenheiten bei der Betrachtung einzubeziehen. Allgemien gilt jedoch, ist die Bebauung insgesamt lockerer, kann die Baulücke größer sein als bei typischerweise eher engerer Bebauung.
* **Empfohlener Wert:** 40 – 100 m.
* **Wirkung:**

Der Algorithmus identifiziert Flächen, die zwischen bestehenden Siedlungs-Clusterflächen liegen und mindestens zu 75 % an diese angrenzen.
Nur wenn der Anteil der Berührung mit Siedlungsflächen deutlich größer ist als jener zur umgebenden freien Landschaft, wird die Fläche als Baulücke klassifiziert.
Diese Flächen werden anschließend als integraler Bestandteil des Innenbereichs behandelt und im Polygon zusammengeführt (vgl. Harig 2024, Kap. 4.5.4, Abb. 4.10).

---

## 4. Hinweise zur Parametrisierung

| Parameter                  | Empfohlener Startwert | Typischer Einfluss                | Sensitivität |
| -------------------------- | --------------------: | --------------------------------- | ------------ |
| `min_overlap_blocks`       |                    18 | Bearbeitungszeit reduzieren       | hoch         |
| `global_footprint_density` |              auto (0) | Abgrenzung kl. Siedlungen         | hoch         |
| `min_area`                 |               56,8 m² | Keine kl. Gebäude am Siedlungsrand| gering       |
| `min_bdg_count`            |                    20 | Mindestanzahl Gebäude Ortsteil    | mittel       |
| `min_patch_size`           |             10 000 m² | Rauschunterdrückung               | mittel       |
| `max_hole_size`            |             10 000 m² | Füllung Außenber. im Innenber.    | mittel       |
| `max_gap_size`             |              4 900 m² | Schließen von Baulücken           | hoch         |

---

## 5. Quellen und Referenzen

* Harig, O. (2024): *Automatisierte Abgrenzung von Innenbereichen einschließlich Ergebnisevaluierung – Grundlage für ein Siedlungsflächenmonitoring.* TU Dresden, Fakultät Umweltwissenschaften.

  * Kapitel 4.4 – Methodik
  * Kapitel 4.5 – Parametrisierung des Verfahrens
  * Anhang A.1.3 – Expertenbefragung
* Bukies, M.; Meyer, G.; Rabe, H. (2009): *Abgrenzung des Innenbereichs im unbeplanten Siedlungsgebiet.*
* Spannowsky u. a. (2020, 2022): *Baugesetzbuch – Kommentar zu § 34 BauGB.*
