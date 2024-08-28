TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{{TITLE1}}</title>
  <link rel="stylesheet" href="https://fonts.googleapis.com/icon?family=Material+Icons">
  <link rel="stylesheet" href="https://code.getmdl.io/1.3.0/material.indigo-pink.min.css">
  <script defer src="https://code.getmdl.io/1.3.0/material.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.125/build/three.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/scatter-gl@0.0.13/lib/scatter-gl.min.js"></script>
  <style>
    body {
      font-family: Arial, Helvetica, sans-serif;
      margin: 0;
    }

    #container {
      width: 100vw;
      height: 100vh;
    }

    #controls {
      position: fixed;
      top: 0;
      left: 0;
      padding: 10px;
    }

    #point-info {
      position: fixed;
      bottom: 10px;
      left: 10px;
      padding: 10px;
      background-color: rgba(255, 255, 255, 0.25);
      border-radius: 5px;
      max-width: 400px;
      max-height: 300px;
    }

    #title {
      position: fixed;
      top: 10px;
      left: 50%;
      transform: translateX(-50%);
      padding: 10px;
      background-color: rgba(255, 255, 255, 0.8);
      border-radius: 5px;
      font-size: 24px;
      font-weight: bold;
    }

    .mdl-js-radio {
      margin-left: 10px;
    }

    .button {
      margin-top: 18px;
    }

    #color-legend {
      position: fixed;
      top: 10px;
      right: 10px;
      padding: 10px;
      background-color: rgba(255, 255, 255, 0.8);
      border-radius: 5px;
      display: none;
    }

    .legend-item {
      display: flex;
      align-items: center;
      margin-bottom: 5px;
    }

    .legend-color {
      width: 20px;
      height: 20px;
      margin-right: 10px;
      border-radius: 50%;
    }
  </style>
</head>
<body>
  <div id="container"></div>
  <div id="title">{{TITLE2}}</div>
  <div id="controls">
    <div class="interactions control">
      <h5>Controls</h5>
      <label class="mdl-radio mdl-js-radio" for="pan-interaction">
        <input type="radio" class="mdl-radio__button" id="pan-interaction" name="interactions" value="pan" checked>
        <span class="mdl-radio__label">Pan</span>
      </label>
      <label class="mdl-radio mdl-js-radio" for="select-interaction">
        <input type="radio" class="mdl-radio__button" id="select-interaction" name="interactions" value="select">
        <span class="mdl-radio__label">Select</span>
      </label>
    </div>
    <div class="button control">
      <button id="toggle-orbit" class="mdl-button mdl-js-button mdl-button--raised">
        Toggle Orbit
      </button>
    </div>
    <hr>
    <div class="color control">
      <h5>Color by</h5>
      <select id="color-select" class="mdl-textfield__input">
        <option value="-1">Default</option>
        <!-- Other options will be dynamically added here -->
      </select>
    </div>
    <div id="point-info">
      <pre id="messages"></pre>
    </div>
  </div>
  <div id="color-legend"></div>
  <script>
    // Embedded CSV data
    const dataCSV = `ID,x,y,z
a,1,2,3
b,2,3,4`;

    const metadataCSV = `ID,ch1,ch2,ch3
a,A,B,C
b,B,C,D`;

async function generatePlotHTML() {
      // Parse CSV data
      const dataPoints = parseDataCSV(dataCSV);
      const [pointMetadata, datasetMetadata] = parseMetadataCSV(metadataCSV);

      // Determine point IDs
      const pointIDs = dataPoints.map(point => point.ID);

      // Create color options
      const colorSelect = document.getElementById('color-select');
      pointMetadata.forEach((meta, i) => {
        const option = document.createElement('option');
        option.value = i.toString();
        option.textContent = meta.name;
        colorSelect.appendChild(option);
      });

      const colorLegend = document.getElementById('color-legend');

      colorSelect.addEventListener('change', function() {
        currentColorIndex = parseInt(this.value);
        if (currentColorIndex === -1) {
          scatterGL.setPointColorer(null);
          colorLegend.style.display = 'none';
        } else {
          scatterGL.setPointColorer(createPointColorer(currentColorIndex));
          updateColorLegend(currentColorIndex);
        }
      });

      function createPointColorer(colorIndex) {
        return function(i, selectedIndices, hoverIndex) {
          var labelIndex = datasetMetadata[i][colorIndex];
          var opaque = false;
          if (opaque) {
            return opaqueColorsByLabel[labelIndex];
          } else {
            if (hoverIndex === i) {
              return "red";
            }
            // If nothing is selected, return the heavy color
            if (selectedIndices.size === 0) {
              return heavyTransparentColorsByLabel[labelIndex];
            }
            // Otherwise, keep the selected points heavy and non-selected light
            else {
              var isSelected = selectedIndices.has(i);
              return isSelected
                ? heavyTransparentColorsByLabel[labelIndex]
                : lightTransparentColorsByLabel[labelIndex];
            }
          }
        };
      }

      function updateColorLegend(colorIndex) {
        const meta = pointMetadata[colorIndex];
        colorLegend.innerHTML = `<h5>${meta.name}</h5>`;

        // Sort the values alphabetically
        const sortedValues = [...meta.values].sort((a, b) => a.localeCompare(b));

        sortedValues.forEach((value) => {
          const index = meta.values.indexOf(value);
          const legendItem = document.createElement('div');
          legendItem.className = 'legend-item';
          legendItem.innerHTML = `
            <div class="legend-color" style="background-color: ${heavyTransparentColorsByLabel[index]}"></div>
            <div>${value}</div>
          `;
          colorLegend.appendChild(legendItem);
        });
        colorLegend.style.display = 'block';
      }

      const pointInfoTip = "Select a point for info";
      let clickedPointMetadata = null;
      let currentColorIndex = -1;

      function getPointMetadata(i) {
        var metadata = {};
        for (var j = 0; j < pointMetadata.length; j++) {
          let propertyName = pointMetadata[j].name;
          let propertyValue = pointMetadata[j].values[datasetMetadata[i][j]];
          metadata[propertyName] = propertyValue;
        }
        var metadata_str = Object.entries(metadata)
          .map(([key, value]) => `${key}: ${value}`)
          .join('\\n');
        return `ID: ${pointIDs[i]}\\n${metadata_str}`;
      }

      const dataset = new ScatterGL.Dataset(dataPoints.map(point => [point.x, point.y, point.z]), datasetMetadata);
      const scatterGL = new ScatterGL(document.getElementById("container"), {
        onClick: function(point) {
          if (point === null) {
            clickedPointMetadata = null;
            setMessage(pointInfoTip);
          } else {
            pointClicked = true;
            clickedPointMetadata = getPointMetadata(parseInt(point));
            setMessage(`Clicked point ${clickedPointMetadata}`);
          }
        },
        onHover: function(point) {
          if (point === null) {
            if (clickedPointMetadata === null) {
              setMessage(pointInfoTip);
            } else {
              setMessage(`Clicked point ${clickedPointMetadata}`);
            }
          } else {
            setMessage(`Hovered point ${getPointMetadata(parseInt(point))}`);
          }
        },
        onSelect: function(points) {
          var message = "";
          if (points.length === 0) {
            message = "no selection";
          } else if (points.length === 1) {
            message = "selected " + points;
          } else {
            message = "selected " + points.length + " points";
          }
          setMessage(message);
        },
        renderMode: "POINT",
        orbitControls: {
          zoomSpeed: 1.125,
        },
      });
      scatterGL.render(dataset);

      window.addEventListener("resize", function() {
        scatterGL.resize();
      });

      document.querySelectorAll('input[name="interactions"]').forEach(function(inputElement) {
        inputElement.addEventListener("change", function() {
          if (inputElement.value === "pan") {
            scatterGL.setPanMode();
          } else if (inputElement.value === "select") {
            scatterGL.setSelectMode();
          }
        });
      });

      var hues = [...Array(10)].map((_, i) => Math.floor((255 / 10) * i));
      var lightTransparentColorsByLabel = hues.map(
        hue => `hsla(${hue}, 100%, 50%, 0.05)`
      );
      var heavyTransparentColorsByLabel = hues.map(
        hue => `hsla(${hue}, 100%, 50%, 0.75)`
      );
      var opaqueColorsByLabel = hues.map(hue => `hsla(${hue}, 100%, 60%, 1)`);

      document.querySelectorAll('input[name="color"]').forEach(function(inputElement) {
        inputElement.addEventListener("change", function() {
          currentColorIndex = parseInt(inputElement.value);
          if (currentColorIndex === -1) {
            scatterGL.setPointColorer(null);
          } else {
            scatterGL.setPointColorer(function(i, selectedIndices, hoverIndex) {
              var labelIndex = datasetMetadata[i][currentColorIndex];
              var opaque = false;
              if (opaque) {
                return opaqueColorsByLabel[labelIndex];
              } else {
                if (hoverIndex === i) {
                  return "red";
                }
                // If nothing is selected, return the heavy color
                if (selectedIndices.size === 0) {
                  return heavyTransparentColorsByLabel[labelIndex];
                }
                // Otherwise, keep the selected points heavy and non-selected light
                else {
                  var isSelected = selectedIndices.has(i);
                  return isSelected
                    ? heavyTransparentColorsByLabel[labelIndex]
                    : lightTransparentColorsByLabel[labelIndex];
                }
              }
            });
          }
        });
      });

      const toggleOrbitButton = document.getElementById('toggle-orbit');
      toggleOrbitButton.addEventListener('click', () => {
        if (scatterGL.isOrbiting()) {
          scatterGL.stopOrbitAnimation();
        } else {
          scatterGL.startOrbitAnimation();
        }
      });

      function setMessage(message) {
        console.log(message);
        document.getElementById("messages").innerHTML = message;
      }
      setMessage(pointInfoTip);
    }

    function parseDataCSV(csv) {
      const rows = csv.trim().split('\\n');
      const headers = rows[0].split(',');

      return rows.slice(1).map(row => {
        const values = row.split(',');
        const point = {};
        headers.forEach((header, index) => {
          point[header] = header === 'ID' ? values[index] : parseFloat(values[index]);
        });
        return point;
      });
    }

    function parseMetadataCSV(csv) {
      const rows = csv.trim().split('\\n');
      const headers = rows[0].split(',');

      const metadataByID = {};
      rows.slice(1).forEach(row => {
        const values = row.split(',');
        metadataByID[values[0]] = {};
        headers.slice(1).forEach((header, index) => {
          metadataByID[values[0]][header] = values[index + 1];
        });
      });

      const pointMetadata = headers.slice(1).map(header => ({
        name: header,
        values: [...new Set(Object.values(metadataByID).map(item => item[header]))]
      }));

      const metadataLookup = pointMetadata.reduce((acc, meta) => {
          acc[meta.name] = meta.values.reduce((acc2, value, i) => {
              acc2[value] = i;
              return acc2;
          }, {});
          return acc;
      }, {});

      const datasetMetadata = Object.values(metadataByID).map(meta => {
        return pointMetadata.map(prop => metadataLookup[prop.name][meta[prop.name]]);
      });

      return [pointMetadata, datasetMetadata];
    }

    generatePlotHTML();
  </script>
</body>
</html>
"""


def generate_scattergl_html(
    title1: str,
    title2: str,
    data_csv_content: str,
    metadata_csv_content: str,
):
    """Generate an interactive 3D scatter plot using scatter-gl."""
    # Replace placeholders in the template
    html_content = TEMPLATE.replace("{{TITLE1}}", title1)
    html_content = html_content.replace("{{TITLE2}}", title2)
    html_content = html_content.replace(
        "ID,x,y,z\\na,1,2,3\\nb,2,3,4", data_csv_content
    )
    html_content = html_content.replace(
        "ID,ch1,ch2,ch3\\na,A,B,C\\nb,B,C,D", metadata_csv_content
    )

    return html_content
