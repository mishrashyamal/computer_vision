# computer_vision


# Doc Scanner using Computer Vision in MATLAB

## Overview

This project implements a document scanner using MATLAB and Computer Vision techniques. The algorithm extracts edges, detects lines using the Hough Transform, identifies intersections, and performs image rectification to transform the document into a properly aligned format.

## Features

- **Edge Detection:** Uses edge detection algorithms to identify prominent edges in the document.
- **Hough Transform:** Detects straight lines in the image to facilitate document rectification.
- **Line Filtering:** Filters and selects the most relevant lines for rectification.
- **Intersection Detection:** Identifies intersection points of the detected lines for geometric transformations.
- **Image Rectification:** Applies a homography transformation to align the document and output a scanned version.

## How It Works

1. **Edge Detection:** The algorithm applies a Gaussian filter for smoothing and then detects edges using gradient operators.
2. **Hough Transform:** It computes the Hough Transform to detect lines in the edge-detected image.
3. **Line Selection:** Relevant lines are selected based on local maxima in the Hough Transform.
4. **Intersection Detection:** Intersection points of the selected lines are computed to define the region for rectification.
5. **Image Rectification:** Using the detected intersection points, the image is warped to produce a top-down view of the document.

## Prerequisites

- MATLAB (with Computer Vision Toolbox)
- Image Processing functions in MATLAB

## Running the Project

1. **Load the Image**  
   The script starts by loading an image (`test.png`). Make sure to replace `'test.png'` with the path to your image.

2. **Edge Detection & Hough Transform**  
   It computes edges and applies the Hough Transform to detect lines in the document.

3. **Line Filtering & Intersection Detection**  
   Relevant lines are selected, and their intersection points are computed.

4. **Rectification**  
   The document is rectified to produce a clean scanned version.

To run the code, execute the `doc_scanner.m` script in MATLAB.

```bash
doc_scanner.m

