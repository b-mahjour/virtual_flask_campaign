import React, { useRef, useEffect } from "react";
import colormap from "colormap";

const HeatmapCanvas = ({ data, setter, sel, classMap }) => {
  const canvasRef = useRef(null);

  useEffect(() => {
    console.log(classMap);
    if (Object.keys(classMap).length === 0) {
      return;
    }
    // const sortedKeys = Object.keys(data).sort();
    let groupedKeys = [];
    if (classMap["A"] || classMap["B"]) {
      ["A", "B"].forEach((className) => {
        if (classMap[className]) {
          const keys = classMap[className].filter((key) => key in data);
          keys.sort(); // Sort keys alphabetically within each class
          groupedKeys.push({ className, keys });
        }
      });
    }
    Object.keys(classMap).forEach((className) => {
      if (className !== "A" && className !== "B") {
        const keys = classMap[className].filter((key) => key in data);
        keys.sort(); // Sort keys alphabetically within each class
        groupedKeys.push({ className, keys });
      }
    });

    const canvas = canvasRef.current;
    const ctx = canvas.getContext("2d");

    const cellWidth = 20;
    const cellHeight = 50;
    const dx = 2;
    const dy = 2;
    const x0 = 25;
    const y0 = 25;
    const groupGap = 10; // Gap between groups on the x-axis

    const numRows = 5; // Since each array has size 5
    const numCols = groupedKeys.reduce(
      (sum, group) => sum + group.keys.length,
      0
    );

    canvas.width = cellWidth * numCols + dx * (numCols + 1) + 100;
    canvas.height = cellHeight * numRows + dy * (numRows + 1) + 50;
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Find the maximum value in the data for normalization
    let maxValue = -Infinity;
    let minValue = Infinity;
    Object.values(data).forEach((values) => {
      values.forEach((value) => {
        if (value > maxValue) {
          maxValue = value;
        }
        if (value < minValue) {
          minValue = value;
        }
      });
    });

    if (maxValue === 0) maxValue = 1;

    const nshades = Math.max(9, maxValue - minValue + 1);
    const colors = colormap({
      colormap: "viridis",
      nshades: nshades,
      format: "hex",
      alpha: 1,
    });

    const getColor = (value) => {
      const normalizedValue = (value - minValue) / (maxValue - minValue);
      const colorIndex = Math.round(normalizedValue * (nshades - 1));
      return colors[colorIndex];
    };

    const cellPositions = [];
    let currentColIndex = 0;
    let currentGroupIndex = 0;

    groupedKeys.forEach(({ className, keys }, groupIndex) => {
      if (groupIndex > 1) {
        currentGroupIndex += 1; // Add one column gap between groups
      }

      // Draw class label
      ctx.fillStyle = "black";
      ctx.font = "8px Arial";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      const classLabelX =
        currentColIndex * cellWidth +
        (currentColIndex + keys.length) * dx +
        x0 +
        (keys.length * cellWidth) / 2 +
        currentGroupIndex * groupGap;
      ctx.fillText(className, classLabelX, y0 - 4);

      // Draw the heatmap cells for each key in the class
      keys.forEach((key) => {
        const values = data[key];
        for (let rowIndex = 0; rowIndex < numRows; rowIndex++) {
          const value = values[rowIndex];

          const x =
            currentColIndex * cellWidth +
            (currentColIndex + 1) * dx +
            x0 +
            currentGroupIndex * groupGap;
          const y = rowIndex * cellHeight + (rowIndex + 1) * dy + y0;

          ctx.fillStyle = getColor(value);
          ctx.fillRect(x, y, cellWidth, cellHeight);

          // Highlight selected cells
          ctx.strokeStyle = "black";
          if (sel[rowIndex] !== null && sel[rowIndex].key === key) {
            ctx.strokeStyle = "red";
          }
          ctx.lineWidth = 1;
          ctx.strokeRect(x, y, cellWidth, cellHeight);

          // Draw value text in cell
          ctx.font = "14px Arial";
          ctx.textAlign = "center";
          ctx.textBaseline = "middle";
          ctx.fillStyle = "black";
          ctx.fillText(value, x + cellWidth / 2, y + cellHeight / 2);

          cellPositions.push({
            x,
            y,
            width: cellWidth,
            height: cellHeight,
            colIndex: currentColIndex,
            rowIndex,
            key,
            value,
          });
        }
        currentColIndex++;
      });
    });

    ctx.fillStyle = "black";
    ctx.font = "14px Arial";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";

    // Draw column labels (keys)
    currentColIndex = 0;
    currentGroupIndex = 0;
    groupedKeys.forEach(({ keys }, groupIndex) => {
      if (groupIndex > 0) {
        currentGroupIndex += 1; // Add one column gap between groups
      }

      keys.forEach((key) => {
        ctx.fillStyle = "black";
        ctx.font = "14px Arial";
        ctx.textAlign = "center";
        ctx.textBaseline = "top";
        ctx.fillText(
          key,
          currentColIndex * cellWidth +
            cellWidth / 2 +
            (currentColIndex + 1) * dx +
            x0 +
            currentGroupIndex * groupGap,
          numRows * cellHeight + 5 + (numRows + 1) * dy + y0
        );
        currentColIndex++;
      });
    });

    // Draw row labels
    ctx.textAlign = "left";
    for (let rowIndex = 0; rowIndex < numRows; rowIndex++) {
      ctx.fillText(
        `d${rowIndex + 1}`,
        5,
        rowIndex * cellHeight + cellHeight / 2 + y0 + (rowIndex + 1) * dy
      );
    }
    // for (let colIndex = 0; colIndex < numCols; colIndex++) {
    //   const key = sortedKeys[colIndex];
    //   ctx.fillText(
    //     key,
    //     colIndex * cellWidth + cellWidth / 2 + (colIndex + 1) * dx + x0,
    //     numRows * cellHeight + 5 + (numRows + 1) * dy + y0
    //   );
    // }

    // // add left hand row labels
    // ctx.textAlign = "left";
    // for (let rowIndex = 0; rowIndex < numRows; rowIndex++) {
    //   ctx.fillText(
    //     `d${rowIndex + 1}`,
    //     7,
    //     rowIndex * cellHeight + cellHeight / 2 - 7 + y0 + (rowIndex + 1) * dy
    //   );
    // }

    const handleClick = (e) => {
      const rect = canvas.getBoundingClientRect();
      const scaleX = canvas.width / rect.width;
      const scaleY = canvas.height / rect.height;
      const x = (e.clientX - rect.left) * scaleX;
      const y = (e.clientY - rect.top) * scaleY;

      const cell = cellPositions.find(
        (cell) =>
          x >= cell.x &&
          x < cell.x + cell.width &&
          y >= cell.y &&
          y < cell.y + cell.height
      );

      if (cell) {
        let new_sel = [...sel];
        if (
          new_sel[cell.rowIndex] !== null &&
          new_sel[cell.rowIndex].key === cell.key
        ) {
          new_sel[cell.rowIndex] = null;
        } else {
          new_sel[cell.rowIndex] = cell;
        }
        setter(new_sel);
      }
    };

    canvas.addEventListener("click", handleClick);

    return () => {
      canvas.removeEventListener("click", handleClick);
    };
  }, [data, setter, sel]);

  return (
    <div style={{ display: "flex", width: "100%" }}>
      <canvas ref={canvasRef} style={{ backgroundColor: "lightgrey" }} />
    </div>
  );
};

export default HeatmapCanvas;
