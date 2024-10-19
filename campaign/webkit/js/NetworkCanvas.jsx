import React, { useState, useEffect, useRef } from "react";
import ReactDOM from "react-dom";
import "../static/main.css";

const NetworkCanvas = ({ data, selectedNode, setSelectedNode }) => {
  const canvasRef = useRef(null);
  const parentRef = useRef(null);
  const [dimensions, setDimensions] = useState({ width: 2500, height: 1000 }); // Default size
  const [offset, setOffset] = useState({ x: 0, y: 0 });
  const [zoom, setZoom] = useState(1);
  const nodeRadius = 4;

  useEffect(() => {
    const handleKeyDown = (event) => {
      const { key } = event;
      switch (key) {
        case "w":
        case "ArrowUp":
          setOffset((prevOffset) => ({ ...prevOffset, y: prevOffset.y + 50 }));
          break;
        case "s":
        case "ArrowDown":
          setOffset((prevOffset) => ({ ...prevOffset, y: prevOffset.y - 50 }));
          break;
        case "a":
        case "ArrowLeft":
          setOffset((prevOffset) => ({ ...prevOffset, x: prevOffset.x + 50 }));
          break;
        case "d":
        case "ArrowRight":
          setOffset((prevOffset) => ({ ...prevOffset, x: prevOffset.x - 50 }));
          break;
        case "z": // Zoom in
          setZoom((prev) => prev * 1.1);
          break;
        case "x": // Zoom out
          setZoom((prev) => prev * 0.9);
          break;
        default:
          break;
      }
    };

    window.addEventListener("keydown", handleKeyDown);
    return () => window.removeEventListener("keydown", handleKeyDown);
  }, []);

  useEffect(() => {
    const observeTarget = parentRef.current;
    const resizeObserver = new ResizeObserver((entries) => {
      if (entries[0].contentRect) {
        setDimensions({
          width: entries[0].contentRect.width,
          height: entries[0].contentRect.height,
        });
      }
    });

    if (observeTarget) {
      resizeObserver.observe(observeTarget);
    }

    // Cleanup
    return () => {
      if (observeTarget) {
        resizeObserver.unobserve(observeTarget);
      }
    };
  }, []);

  useEffect(() => {
    const drawCanvas = () => {
      if (!data || data === undefined) {
        return;
      }
      if (!data.links) {
        return;
      }
      const canvas = canvasRef.current;
      const context = canvas.getContext("2d");

      context.clearRect(0, 0, dimensions.width, dimensions.height);
      data.links.forEach((link) => {
        const sourceNode = data.nodes.find(
          (node) => node.sm === link.source_node_smiles
        );
        const targetNode = data.nodes.find(
          (node) => node.sm === link.destination_node_smiles
        );
        if (sourceNode && targetNode) {
          context.beginPath();
          context.moveTo(
            sourceNode.x * dimensions.width + offset.x,
            sourceNode.y * dimensions.height + offset.y
          );
          context.lineTo(
            targetNode.x * dimensions.width + offset.x,
            targetNode.y * dimensions.height + offset.y
          );
          context.strokeStyle = "rgba(0, 0, 0, 1)";
          context.lineWidth = 1;
          context.globalAlpha = 1;
          const link_array = [
            link.source_node_smiles,
            link.destination_node_smiles,
          ];
          if (selectedNode === null || selectedNode === undefined) {
            context.globalAlpha = 1;
          } else {
            if (link_array.includes(selectedNode.sm)) {
              context.globalAlpha = 1;
            } else {
              context.globalAlpha = 0.15;
            }
          }
          context.stroke();
        }
      });

      data.nodes.forEach((node) => {
        const x = node.x * dimensions.width + offset.x;
        const y = node.y * dimensions.height + offset.y;

        // const scale = 0.2;
        // const img = new Image();

        context.beginPath();
        context.moveTo(x + nodeRadius, y + nodeRadius);
        context.arc(x, y, nodeRadius, 0, 2 * Math.PI);
        context.lineWidth = 1;

        context.fillStyle = "#f5f5f5";
        if (node.product === true) {
          // console.log(node.product, node)
          context.strokeStyle = "#FFC43D";
        } else {
          context.strokeStyle = "#3d0814";
        }

        if (node.reaction_hits.length > 0 || node.propagations == 0) {
          context.fillStyle = "#748067";
        }
        // console.log(node.target_found)
        if (node.target_found === true) {
          // console.log("hello?")
          context.fillStyle = "#85F942";
        }

        if (selectedNode === null || selectedNode === undefined) {
          context.globalAlpha = 1;
        } else {
          if (selectedNode.id === node.id) {
            context.strokeStyle = "#FF0000";
            context.globalAlpha = 1;
          } else {
            context.globalAlpha = 0.15;
          }
        }

        context.fill();
        context.stroke();

        // img.onload = () => {
        //   const scaledWidth = img.naturalWidth * scale * zoom;
        //   const scaledHeight = img.naturalHeight * scale * zoom;
        //   const adjustedX = x - scaledWidth / 2;
        //   const adjustedY = y - scaledHeight / 2;
        //   context.globalAlpha = 1;
        //   if (selectedNode === null || selectedNode === undefined) {
        //     context.globalAlpha = 1;
        //   } else{
        //     if (selectedNode.id === node.id) {
        //       context.globalAlpha = 1;
        //     } else {
        //       context.globalAlpha = 0.15;
        //     }
        //   }
        //   context.drawImage(
        //     img,
        //     adjustedX,
        //     adjustedY,
        //     scaledWidth,
        //     scaledHeight
        //   );
        //   // console.log(adjustedX, adjustedY, scaledWidth, scaledHeight);
        //   context.beginPath();

        //   node.hull.forEach((point, index) => {
        //     const hullX = adjustedX + point[0] * scale * zoom;
        //     const hullY = adjustedY + point[1] * scale * zoom;
        //     if (index === 0) {
        //       context.moveTo(hullX, hullY);
        //     } else {
        //       context.lineTo(hullX, hullY);
        //     }
        //   });

        //   //         // Connect back to the first point to close the polygon
        //   if (node.hull.length > 0) {
        //     const firstHullPoint = node.hull[0];
        //     const hullX = adjustedX + firstHullPoint[0] * scale * zoom;
        //     const hullY = adjustedY + firstHullPoint[1] * scale * zoom;
        //     context.lineTo(hullX, hullY);
        //   }

        //   if (node.product) {
        //     context.strokeStyle = "#1f3569";
        //   } else {
        //     context.strokeStyle = "#FF0000";
        //   }

        //   // if (selectedNode && node.id === selectedNode.id && node.hull) {
        //   //   context.strokeStyle = "#FF0000"; // Example: Red border
        //   // } else {
        //   //   context.strokeStyle = "#000000"; // Example: Black border
        //   //   // context.strokeStyle = interpolateViridis(10 / 10);
        //   // }
        //   context.lineWidth = 2; // Example: Border width
        //   context.stroke();
        //   context.closePath();
        // };
        // // };
        // img.src = `data:image/png;base64,${node.img}`;
      });
    };

    drawCanvas();
  }, [data, dimensions, offset, zoom, selectedNode]);

  useEffect(() => {
    const canvas = canvasRef.current;
    const handleClick = (event) => {
      if (!data) {
        return;
      }
      if (!data.links) {
        return;
      }
      // console.log(selectedNode);
      if (selectedNode !== null) {
        setSelectedNode(null);
        return;
      }
      const rect = canvas.getBoundingClientRect();

      const scaleX = canvas.width / rect.width; // Relationship bitmap vs. element for X
      const scaleY = canvas.height / rect.height; // Relationship bitmap vs. element for Y

      const canvasX = (event.clientX - rect.left) * scaleX;
      const canvasY = (event.clientY - rect.top) * scaleY;

      // Adjust click coordinates based on the zoom and offset
      const adjustedX = canvasX - offset.x;
      const adjustedY = canvasY - offset.y;

      // Find if a node was clicked
      const clickedNode = data.nodes.find((node) => {
        const nodeX = node.x * dimensions.width;
        const nodeY = node.y * dimensions.height;
        return (
          Math.sqrt((adjustedX - nodeX) ** 2 + (adjustedY - nodeY) ** 2) <
          nodeRadius
        );
      });
      // console.log("selected node", clickedNode);
      // console.log(setSelectedNode);
      setSelectedNode(clickedNode);
    };

    canvas.addEventListener("click", handleClick);

    return () => canvas.removeEventListener("click", handleClick);
  }, [dimensions, data, offset, zoom]);

  return (
    <div ref={parentRef} style={{ width: "100%", height: "100%" }}>
      <canvas
        ref={canvasRef}
        width={dimensions.width}
        height={dimensions.height}
      />
    </div>
  );
};

export default NetworkCanvas;
