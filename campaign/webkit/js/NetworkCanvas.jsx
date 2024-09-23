import React, { useState, useEffect, useRef } from "react";
import ReactDOM from "react-dom";
import "../static/main.css";

function NetworkCanvas(data) {
  const canvasRef = useRef(null);
  const parentRef = useRef(null);
  const [dimensions, setDimensions] = useState({ width: 2000, height: 1000 }); // Default size
  const [offset, setOffset] = useState({ x: 0, y: 0 });
  const [zoom, setZoom] = useState(1); // Zoom level state

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
      data = data.data;
      console.log(data.links);
      if (!data.links) {
        return;
      }
      // if (canvasRef.length === undefined) {
      //   return;
      // }
      console.log("hi", parentRef)
      const canvas = canvasRef.current;
      const context = canvas.getContext("2d");

      context.clearRect(0, 0, dimensions.width, dimensions.height);
      console.log("hello1!");
      data.links.forEach((link) => {
        const sourceNode = data.nodes.find((node) => node.id === link[0]);
        const targetNode = data.nodes.find((node) => node.id === link[1]);
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
          context.strokeStyle = "rgba(0, 0, 0, 0.25)"; // Semi-transparent lines
          context.lineWidth = 2; // Set line width

          const link_array = [link[0], link[1]];
          // if (selectedNode) {
          //   // console.log(link_array, selectedNode)
          //   if (link_array.includes(selectedNode.id)) {
          //     context.globalAlpha = 1;
          //   } else {
          //     context.globalAlpha = 0.15;
          //   }
          // } else {
          //   context.globalAlpha = 1;
          // }
          context.globalAlpha = 1;

          context.stroke();
        }
      });

      // function getConnectedNodeIds(selectedNodeId, links) {
      //   return links.reduce((connectedNodes, link) => {
      //     if (link.source === selectedNodeId) {
      //       connectedNodes.push(link.target);
      //     } else if (link.target === selectedNodeId) {
      //       connectedNodes.push(link.source);
      //     }
      //     return connectedNodes;
      //   }, []);
      // }
      console.log("hello2!");

      data.nodes.forEach((node) => {
        // const connectedNodeIds = selectedNode
        //   ? getConnectedNodeIds(selectedNode.id, data.links)
        //   : [];
        // const isNodeSelectedOrConnected =
        //   node.id === selectedNode?.id || connectedNodeIds.includes(node.id);
        // console.log(dimensions);
        const x = node.x * dimensions.width + offset.x;
        const y = node.y * dimensions.height + offset.y;

        const img = new Image();
        img.onload = () => {
          const scaledWidth = img.naturalWidth * 0.6 * zoom;
          const scaledHeight = img.naturalHeight * 0.6 * zoom;
          const adjustedX = x - scaledWidth / 2;
          const adjustedY = y - scaledHeight / 2;
          // if (selectedNode) {
          //   context.globalAlpha = isNodeSelectedOrConnected ? 1 : 0.15;
          // } else {
          //   context.globalAlpha = 1;
          // }
          context.globalAlpha = 1;
          context.drawImage(
            img,
            adjustedX,
            adjustedY,
            scaledWidth,
            scaledHeight
          );
          console.log(adjustedX, adjustedY, scaledWidth, scaledHeight);
          context.beginPath();

          node.hull.forEach((point, index) => {
            const hullX = adjustedX + point[0] * 0.6 * zoom;
            const hullY = adjustedY + point[1] * 0.6 * zoom;
            if (index === 0) {
              context.moveTo(hullX, hullY);
            } else {
              context.lineTo(hullX, hullY);
            }
          });

          //         // Connect back to the first point to close the polygon
          if (node.hull.length > 0) {
            const firstHullPoint = node.hull[0];
            const hullX = adjustedX + firstHullPoint[0] * 0.6 * zoom;
            const hullY = adjustedY + firstHullPoint[1] * 0.6 * zoom;
            context.lineTo(hullX, hullY);
          }

          if (node.product) {
            context.strokeStyle = "#1f3569"; 

          } else{
            context.strokeStyle = "#FF0000";

          }

          // if (selectedNode && node.id === selectedNode.id && node.hull) {
          //   context.strokeStyle = "#FF0000"; // Example: Red border
          // } else {
          //   // context.strokeStyle = "#000000"; // Example: Black border
          //   context.strokeStyle = interpolateViridis(10 / 10);
          // }
          context.lineWidth = 2; // Example: Border width
          context.stroke();
          context.closePath();
        };
        // };
        img.src = `data:image/png;base64,${node.img}`;
      });
    };

    drawCanvas();
  }, [data, dimensions, offset, zoom]);

  return (
    <div ref={parentRef} style={{ width: "100%", height: "100%" }}>
      <canvas
        ref={canvasRef}
        width={dimensions.width}
        height={dimensions.height}
      />
    </div>
  );
}

export default NetworkCanvas;
