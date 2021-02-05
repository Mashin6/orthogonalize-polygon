# Orthogonalize Polygon in python
This script performs squaring/orthogonalization of polygons, in other words making all its angles 90˚ or 180˚, which is usefull for automated adjustment of building footprint shapes. Script uses GeoPandas and Shapely packages and is an improved implementation of JOSM Orthogonalize function. Supports Polygons with internal holes as well as Multipolygons and a typical polygon takes on Core i7-6660U 2.4GHz about ~0.34s processing time. 
 
**Features**
* Supports Polygons with internal rings/holes.
* Supports MultiPolygons
* Skew angle tolerance to preserve angled features such as bay windows.
* Adjustable range of angles that determines whether polygon edge follows the same direction as previous segment or rather turns 90˚.

**Preserving angled features**

<img src="https://github.com/Mashin6/orthogonalize-polygon/blob/master/Orthogonalize_comparison1.png" width="500">

**Comparison to JOSM**

 <img src="https://github.com/Mashin6/orthogonalize-polygon/blob/master/Orthogonalize_comparison2.png" width="500">
 
 **Example:**
 
<img src="https://github.com/Mashin6/orthogonalize-polygon/raw/master/Orthogonalized.png" width="500">
