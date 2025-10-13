# Dokumentation des Tools zur automatisierten Siedlungsabgrenzung

*(nach Harig 2024, TU Dresden – Automatisierte Abgrenzung von Innenbereichen einschließlich Ergebnisevaluierung)*

---

## 1. Einleitung

Das vorliegende Tool dient der **automatisierten Abgrenzung von Siedlungs- bzw. Innenbereichen** auf Basis geometrischer und topographischer Objektdaten.
Es wurde im Rahmen der Dissertation von **Oliver Harig (2024)** entwickelt und zielt darauf ab, **einheitliche, datenbasierte und nachvollziehbare Siedlungsgrenzen** zu erzeugen, die den planerischen Innenbereich nach § 34 BauGB approximieren.

Die Methodik beruht auf der Kombination von:

* **Gebäudegrundrissen (HU-DE)** als Indikator realer Bebauung,
* **topographischen Linien (ATKIS-Straßennetz)** zur Partitionierung der Raumstruktur,
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

  * Bei **hohen Werten** (z. B. > 0.25) werden nur kompakte, geschlossene Baublöcke als Siedlung erkannt.
  * Bei **niedrigen Werten** (z. B. < 0.15) werden auch lockere Bebauungen (z. B. Dorfkerne) erfasst.
* **Empfohlener Bereich:** 0.18 – 0.22 (vgl. Kap. 4.5.3, empirische Kalibrierung).
* **Einfluss:** Bestimmt direkt die Mindestdichte für die Aufnahme in den Siedlungszusammenhang und beeinflusst somit die resultierende Fläche der Innenbereiche.

---

### 3.2 **`global_footprint_density` – Globale Gebäudedichte**

* **Typ:** Float
* **Funktion:** Wird als Referenzschwellwert für die Gesamtfläche verwendet, wenn keine lokalen Dichtewerte vorliegen (`calc_footprint_density()`).
* **Definition:** Mittelwert der Gebäudeflächendichte aller Partitionen eines Untersuchungsgebiets.
* **Hintergrund:** Dient der Robustheit der Analyse in Gebieten mit geringer Bebauungsdichte oder unvollständigen Gebäudedaten (vgl. Harig 2024, Kap. 4.5.3).
* **Wirkung:**

  * Niedrige Werte (< 0.1) → höhere Sensitivität, Erfassung auch locker bebauter Gebiete.
  * Hohe Werte (> 0.2) → nur verdichtete Siedlungskerne werden erkannt.
* **Empfohlene Einstellung:** Automatische Berechnung (0 = auto); manuell nur bei stark heterogenen Untersuchungsgebieten setzen.
* **Bemerkung:** Wird intern als relativer UG-Schwellwert zur lokalen Dichte normalisiert.

---

### 3.3 **`min_area` – Mindestgebäudefläche**

* **Typ:** Float
* **Funktion:** Filtert in `input_hu_filter()` zu kleine Gebäudeobjekte aus.
* **Definition:** Gebäude mit einer Grundfläche unterhalb dieses Werts werden bei der Dichteberechnung ignoriert.
* **Ziel:** Ausschluss nicht siedlungsrelevanter Strukturen (z. B. Gartenhäuser, Carports).
* **Hintergrund:** In der Dissertation wird dieser Schritt als „semantische Filterung“ bezeichnet (Kap. 4.4.4 und 4.5.2). Grundlage bildet die Unterscheidung nach Gebäudefunktion (Wohn-, Misch-, Gewerbebau etc.).
* **Empfohlene Werte:** 20 – 30 m², je nach Datenquelle (HU-DE, OSM).
* **Wirkung:** Verringert Rauschen und erhöht die Aussagekraft der Dichteanalyse.

---

### 3.4 **`min_bdg_count` – Minimale Gebäudeanzahl**

* **Typ:** Integer
* **Funktion:** Definiert die Mindestanzahl von Gebäuden, ab der ein Cluster als gültiger Siedlungsbereich gewertet wird.
* **Verwendung:** Sowohl bei der Aggregation (Clusterbildung) als auch in der Post-Filterung.
* **Hintergrund:** Laut empirischer Auswertung (Harig 2024, Kap. 4.5.6) und Expertenbefragung (Anhang A.1.3) gilt ein Ortsteil erst ab etwa **20 Gebäuden oder 1 ha** Fläche als eigenständiger Innenbereich.
* **Empfohlene Werte:** 15 – 25 Gebäude.
* **Wirkung:**

  * Verhindert Klassifikation isolierter Gebäudegruppen als Siedlung.
  * Stabilisiert Ergebnisse in ländlichen Räumen.

---

### 3.5 **`min_patch_size` – Minimale Patch-Größe**

* **Typ:** Float
* **Funktion:** Wird in `patch_remove()` zur Entfernung kleiner Siedlungsinseln eingesetzt.
* **Definition:** Unterhalb dieses Flächenwerts wird eine zusammenhängende Siedlung als Rauschen verworfen.
* **Hintergrund:** Entspricht den empirisch ermittelten Mindestgrößen eines Ortsteils (Kap. 4.5.6).
* **Empfohlener Wert:** 10 000 m² (= 1 ha).
* **Wirkung:**

  * Kleine Fragmente werden entfernt → klarere Innenbereichsgrenzen.
  * Reduziert Übersegmentierung bei dünn besiedelten Regionen.

---

### 3.6 **`max_hole_size` – Maximale Lochgröße**

* **Typ:** Float
* **Funktion:** In `gap_close()` zur Füllung unbebauter Innenflächen innerhalb zusammenhängender Siedlungen.
* **Definition:** Gibt an, bis zu welcher Fläche Freiräume innerhalb eines Siedlungspolygons noch als *Baulücken* gelten und gefüllt werden dürfen.
* **Hintergrund:** Entspricht planerischen Definitionen des Innenbereichs; Baulücken ≤ 4900 m² (≈ 70 × 70 m) gelten als Bestandteil des Innenbereichs (vgl. Harig 2024, Kap. 4.5.4, S. 41 f.).
* **Empfohlener Wert:** 5000 m².
* **Wirkung:**

  * Kleinere Freiflächen (Gärten, Innenhöfe) werden integriert.
  * Verhindert übermäßige Fragmentierung der Siedlungsflächen.

---

### 3.7 **`max_gap_size` – Maximale Lückengröße**

* **Typ:** Float
* **Funktion:** In `gap_close()` zur Verbindung benachbarter Siedlungscluster.
* **Definition:** Bestimmt den maximalen Abstand, bis zu dem zwei Siedlungsflächen verbunden werden dürfen.
* **Hintergrund:** Laut empirischer Erhebung unter Fachplanern liegt die als „Bebauungslücke“ wahrgenommene Grenze bei etwa **70 m** (arithm. Mittel 71,8 m, Anhang A.1.3, Frage 8).
* **Empfohlener Wert:** 60 – 80 m.
* **Wirkung:**

  * Schließt kleinere Abstände zwischen Straßenfronten oder Gebäudereihen.
  * Fördert zusammenhängende Siedlungsflächen in verdichteten Gebieten.

---

## 4. Hinweise zur Parametrisierung

| Parameter                  | Empfohlener Startwert | Typischer Einfluss              | Sensitivität |
| -------------------------- | --------------------: | ------------------------------- | ------------ |
| `min_overlap_blocks`       |                  0.20 | Dichteempfindlichkeit           | hoch         |
| `global_footprint_density` |              auto (0) | globale Normierung              | mittel       |
| `min_area`                 |                 25 m² | Filterwirkung Gebäude           | gering       |
| `min_bdg_count`            |                    20 | Mindestgröße Cluster            | hoch         |
| `min_patch_size`           |             10 000 m² | Rauschunterdrückung             | mittel       |
| `max_hole_size`            |              5 000 m² | Füllung Innenräume              | mittel       |
| `max_gap_size`             |                  70 m | Verbindung angrenzender Flächen | hoch         |

### Wechselwirkungen

* **`min_overlap_blocks`** und **`global_footprint_density`** bestimmen die Dichteempfindlichkeit – eine Anpassung des einen Parameters erfordert meist eine Gegenkorrektur des anderen.
* **`max_hole_size`** und **`max_gap_size`** wirken gemeinsam auf die räumliche Kohärenz der Polygone: zu hohe Werte führen zu Überverschmelzung, zu niedrige zu Fragmentierung.
* **`min_bdg_count`** und **`min_patch_size`** definieren die Mindestgröße von Ortsteilen und sind insbesondere in ländlichen Räumen kritisch für realistische Ergebnisse.

---

## 5. Quellen und Referenzen

* Harig, O. (2024): *Automatisierte Abgrenzung von Innenbereichen einschließlich Ergebnisevaluierung – Grundlage für ein Siedlungsflächenmonitoring.* TU Dresden, Fakultät Umweltwissenschaften.

  * Kapitel 4.4 – Methodik
  * Kapitel 4.5 – Parametrisierung des Verfahrens
  * Anhang A.1.3 – Expertenbefragung
* Bukies, M.; Meyer, G.; Rabe, H. (2009): *Abgrenzung des Innenbereichs im unbeplanten Siedlungsgebiet.*
* Spannowsky u. a. (2020, 2022): *Baugesetzbuch – Kommentar zu § 34 BauGB.*
