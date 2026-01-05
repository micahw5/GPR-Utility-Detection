
# GPR Utility Detection Tool for ArcGIS Pro

![Version](https://img.shields.io/badge/version-9.2-blue)
![Python](https://img.shields.io/badge/python-3.7+-green)
![License](https://img.shields.io/badge/license-MIT-orange)
![ArcGIS](https://img.shields.io/badge/ArcGIS%20Pro-3.0+-purple)

Automated Ground Penetrating Radar (GPR) utility detection toolbox for ArcGIS Pro. Uses hyperbola detection algorithms to identify buried utilities from GSSI GPR data and automatically georeferences results.

![Sample Detection](docs/sample_detection.png)
*Example: Automated utility detection with confidence scores and GPS integration*

## 🎯 Features

- **Automated Hyperbola Detection** - Uses edge detection and curve fitting to identify buried utilities
- **GPS Integration** - Automatically georeferences detections from .DZG files
- **Quality Metrics** - Provides confidence scores and R² values for each detection
- **Batch Processing** - Process multiple survey files simultaneously
- **Multi-View Radargrams** - Generates enhanced, contrast, edge, and detection visualization panels
- **ArcGIS Integration** - Creates point and polyline shapefiles ready for GIS analysis
- **Depth & Diameter Estimation** - Calculates approximate utility depth and size

## 📋 Requirements

### Software
- ArcGIS Pro 3.0 or later
- Python 3.7+

### Python Libraries
- **Required**: `numpy`, `arcpy` (included with ArcGIS Pro)
- **Recommended**: `scipy`, `matplotlib` (for advanced detection and visualization)

### Data Format
- GSSI DZT/DZG file format (SIR-3000, SIR-4000, StructureScan, etc.)

## 🚀 Installation

### Option 1: Direct Download
1. Download `toolboxBeta_version.pyt` from this repository
2. Open ArcGIS Pro
3. In the Catalog pane, right-click **Toolboxes** → **Add Toolbox**
4. Navigate to and select `toolboxBeta_version.pyt`

### Option 2: Clone Repository
```bash
git clone https://github.com/micahw5/toolboxBeta_version.git
cd toolboxBeta_version
```

### Install Optional Dependencies
```bash
# Using ArcGIS Pro Python environment
conda activate arcgispro-py3
conda install scipy matplotlib
```

## 📖 Usage

### Basic Workflow

1. **Prepare Data**
   - Place all .DZT files in a single folder
   - Ensure corresponding .DZG GPS files are in the same folder (optional but recommended)

2. **Run Tool**
   - Open the toolbox in ArcGIS Pro
   - Run "GPR Utility Detection Tool"
   - Select input folder containing .DZT files
   - Choose output folder for results
   - Adjust parameters as needed

3. **Review Results**
   - Examine generated radargrams in `Radargrams/` subfolder
   - Review detected utilities in attribute table
   - Filter by confidence score (>0.7 recommended for high-quality detections)

### Parameters

| Parameter | Description | Default | Notes |
|-----------|-------------|---------|-------|
| **Input Folder** | Folder containing .DZT files | Required | Also looks for .DZG GPS files |
| **Output Folder** | Where to save results | Required | Creates shapefiles and radargrams |
| **Survey Name** | Prefix for output files | "GPR_Survey" | Used in shapefile names |
| **Velocity** | Subsurface velocity (m/ns) | Auto (0.09) | 0.09 typical for urban soil |
| **Min Confidence** | Detection threshold | 0.58 | Lower = more detections, more false positives |
| **Auto-add to map** | Add results to current map | True | Automatically loads shapefiles |

### Output Files

```
Output_Folder/
├── [SurveyName]_Utilities.shp      # Point features for each utility
├── [SurveyName]_Path.shp           # Survey path polyline
└── Radargrams/
    ├── FILE001_V9_2.png            # 4-panel visualization
    ├── FILE002_V9_2.png
    └── ...
```

### Attribute Table Fields

| Field | Description | Type |
|-------|-------------|------|
| **UtilID** | Unique identifier | Text |
| **Depth_m** | Estimated depth in meters | Double |
| **Diam_m** | Estimated diameter in meters | Double |
| **Confidence** | Detection confidence (0-1) | Double |
| **R_Squared** | Hyperbola fit quality (0-1) | Double |
| **Velocity** | Velocity used for depth calc | Double |
| **SourceFile** | Origin .DZT filename | Text |

## 🎨 Interpreting Results

### Confidence Scores
- **> 0.8** (Red markers): High confidence - strong hyperbola fit
- **0.6 - 0.8** (Orange markers): Medium confidence - good fit with some uncertainty
- **< 0.6** (Yellow markers): Lower confidence - requires careful review

### Best Practices
- ⚠️ **Always review automated detections** - Professional verification recommended
- Use confidence scores to prioritize field investigation
- Cross-reference with utility records and as-builts
- Consider soil conditions and depth limitations
- Validate with potholing before excavation

### Typical Detection Rates
Performance varies based on:
- Utility material (metallic > plastic)
- Depth (shallow > deep)
- Soil type (sandy > clay)
- Data quality and GPR settings

## 🔧 Advanced Configuration

### Custom Velocity Values

Different materials have different velocities:

| Material | Velocity (m/ns) |
|----------|-----------------|
| Air | 0.30 |
| Dry sand | 0.15 |
| Saturated sand | 0.06 |
| Dry clay | 0.12 |
| Saturated clay | 0.06 |
| Concrete | 0.10 |
| Asphalt | 0.08 |

### Troubleshooting

**No detections found:**
- Lower minimum confidence (try 0.50)
- Check velocity setting
- Verify .DZT file quality and format

**Too many false positives:**
- Increase minimum confidence (try 0.65)
- Review radargrams for data quality issues
- Check for surface clutter/noise

**GPS coordinates incorrect:**
- Verify .DZG file is in same folder as .DZT
- Check GPS data quality in DZG file
- Ensure proper NMEA formatting

## 🤝 Contributing

Contributions are welcome! Areas for improvement:

- [ ] Support for additional GPR file formats (Sensors & Software, IDS, Mala)
- [ ] 3D visualization capabilities
- [ ] Machine learning classification (utility type detection)
- [ ] Integration with utility asset management systems
- [ ] Performance optimization for large datasets
- [ ] Ground truth validation tools

### How to Contribute
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📊 Validation & Accuracy

This tool is under active development and validation. Users should:
- Conduct independent validation for their specific applications
- Not rely solely on automated detection for critical excavation decisions
- Follow local utility location regulations and best practices

We welcome validation studies and accuracy reports from the community.

## 📝 Citation

If you use this tool in research or publications, please cite:

```
[Micah Williams/Joe Risk]. (2025). GPR Utility Detection Tool for ArcGIS Pro (Version 9.2) 
[Software]. GitHub. https://github.com/micahw5/gpr-utility-detection
```

## 🔗 Related Resources

- [GSSI File Formats Documentation](https://www.geophysical.com/)
- [ArcGIS Pro Python Reference](https://pro.arcgis.com/en/pro-app/latest/arcpy/)
- [GPR Principles and Applications](https://wiki.seg.org/wiki/Ground_penetrating_radar)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## ⚠️ Disclaimer

This tool is provided "as is" without warranty of any kind. Users are responsible for:
- Validating results before excavation
- Following local regulations for utility location
- Obtaining proper clearances and permits
- Professional review of automated detections

Incorrect utility location can result in serious injury, death, or property damage. Always use multiple verification methods.

## 👤 Authors

**[Micah Williams]**
- GitHub: [micahw5](https://github.com/micahw5)
- LinkedIn: [Your LinkedIn](https://linkedin.com/in/yourprofile)
- Email: mwilliams@nelsonintelligencesolutions.com

 
 **[Joe Risk]**
- GitHub: [micahw5](https://github.com/yourusername)
- LinkedIn: [Your LinkedIn](https://linkedin.com/in/yourprofile)
- Email: jrisk@nelsonintelligencesolutions.com

## 🙏 Acknowledgments

- Built with ArcPy and the ArcGIS Pro Python ecosystem
- Uses scipy for advanced numerical processing
- Inspired by the need for accessible GPR interpretation tools
- Thanks to the open-source GIS and geophysics communities

---

**⭐ If you find this tool useful, please star the repository!**

**🐛 Found a bug?** [Open an issue](https://github.com/micahw5/gpr-utility-detection/issues)

**💡 Have a suggestion?** [Start a discussion](https://github.com/micahw5/gpr-utility-detection/discussions)
