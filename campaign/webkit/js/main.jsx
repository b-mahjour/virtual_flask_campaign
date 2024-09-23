import React, { useState, useEffect, useRef } from "react";
import ReactDOM from "react-dom";
import "../static/main.css";
import NetworkCanvas from "./NetworkCanvas";

function Main() {
  const [network, setNetwork] = useState({ nodes: [], links: [] });

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

  function renderNetwork() {
    return network.nodes.map((node) => {
      return (
        <div className="node">
          <div className="node-name">{node.sm}</div>
          <img
            src={"data:image/png;base64," + node.img}
            alt=""
            style={{
              maxWidth: "100%",
              height: "auto",
              display: "block",
              width: "100px",
            }}
          />
        </div>
      );
    });
  }

  return (
    <div className="container">
      <div className="netHolder">
      <NetworkCanvas data={network} />

      </div>
      {/* <div className="main">{renderNetwork()}</div> */}
    </div>
  );
}

// ========================================

ReactDOM.render(<Main />, document.getElementById("reactEntry"));
