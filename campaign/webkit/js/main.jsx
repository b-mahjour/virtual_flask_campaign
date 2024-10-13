import React, { useState, useEffect, useRef } from "react";
import ReactDOM from "react-dom";
import "../static/main.css";
import NetworkCanvas from "./NetworkCanvas";

function Main() {
  const [network, setNetwork] = useState({ nodes: [], links: [] });
  const [selectedNode, setSelectedNode] = useState(null);

  function getHits() {
    fetch("/api/get_network", {
      method: "GET",
    })
      .then((response) => response.json())
      .then((data) => {
        setNetwork(data.network);
      });
  }

  useEffect(() => {
    getHits();
  }, []);

  function renderSelectedInfo() {
    if (!selectedNode) {
      return null;
    }

    let arr = [];
    for (let i in selectedNode.reaction_hits){
      arr.push(<div>
        <div>{selectedNode.reaction_hits[i].wells}</div>
        <div>{selectedNode.reaction_hits[i].exact_mass} ({selectedNode.reaction_hits[i].mass_found})</div>
      </div>)
    }

    return (
      <div className="selected-node">
        <div className="node-name">{selectedNode.sm}</div>
        <img
          src={"data:image/png;base64," + selectedNode.img}
          alt=""
          style={{
            maxWidth: "100%",
            height: "auto",
            display: "block",
            // width: "100px",
          }}
        />
        <div>{arr}</div>
      </div>
    );
  }

  function setSelectedNodeRandomly() {
    setSelectedNode(
      network.nodes[Math.floor(Math.random() * network.nodes.length)]
    );
  }

  function renderCommands() {
    return (
      <div className="commands">
        <div className="buttonHolder">
          <button onClick={setSelectedNodeRandomly}>Target Node</button>
          <button>Download Array</button>
        </div>
      </div>
    );
  }

  return (
    <div className="container">
      <div className="netHolder">
        <NetworkCanvas
          data={network}
          selectedNode={selectedNode}
          setSelectedNode={setSelectedNode}
        />
      </div>
      <div className="rightSide">
        {renderCommands()}
        {renderSelectedInfo()}
      </div>
    </div>
  );
}

// ========================================

ReactDOM.render(<Main />, document.getElementById("reactEntry"));
