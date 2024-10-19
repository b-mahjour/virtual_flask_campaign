import React, { useState, useEffect, useRef } from "react";
import ReactDOM from "react-dom";
import "../static/main.css";
import NetworkCanvas from "./NetworkCanvas";
import HeatmapCanvas from "./HeatmapCanvas";

function Main() {
  const [network, setNetwork] = useState({ nodes: [], links: [] });
  const [seqReactionsDict, setSeqReactionsDict] = useState({});
  const [invMap, setInvMap] = useState({});
  const [classMap, setClassMap] = useState({});
  const [wellToSequenceCounts, setWellToSequenceCounts] = useState({});
  const [selectedNode, setSelectedNode] = useState(null);
  const [selectedSequenceItems, setSelectedSequenceItems] = useState([
    null,
    null,
    null,
    null,
    null,
  ]);

  function getHits() {
    fetch("/api/get_network", {
      method: "GET",
    })
      .then((response) => response.json())
      .then((data) => {
        setNetwork(data.network);
      });
  }

  function getReactions() {
    fetch("/api/get_reactions", {
      method: "GET",
    })
      .then((response) => response.json())
      .then((data) => {
        setClassMap(data.class_map);
        setInvMap(data.inv_map);
        setSeqReactionsDict(data.action_sequence_to_data);
        setWellToSequenceCounts(data.well_to_sequence_counts);
      });
  }

  useEffect(() => {
    getHits();
  }, []);

  useEffect(() => {
    getReactions();
  }, []);

  function renderNodeSelectedInfo() {
    if (!selectedNode) {
      return <div className="selected-node"></div>;
    }

    let arr = [
      <div className="node-hit-item">
        <div className="node-hit-col">rxn-well</div>
        <div className="node-hit-col">E.M. (found)</div>
      </div>,
    ];
    for (let i in selectedNode.reaction_hits) {
      arr.push(
        <div className="node-hit-item">
          <div className="node-hit-col left">
            {selectedNode.reaction_hits[i].wells}
          </div>
          <div className="node-hit-col">
            {selectedNode.reaction_hits[i].exact_mass} (
            {selectedNode.reaction_hits[i].mass_found})
          </div>
        </div>
      );
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
        <div className="node-hit-data">{arr}</div>
      </div>
    );
  }

  function renderSelectedSequenceItems() {
    let arr = [];
    for (let i in selectedSequenceItems) {
      if (selectedSequenceItems[i] === null) {
        arr.push(<div className="reagentRow">-</div>);
      } else {
        let reag = invMap[selectedSequenceItems[i].key];
        arr.push(
          <div className="reagentRow">
            <div className="reagentItem">{reag.description}</div>
            <div className="reagentItem">
              <img
                src={"data:image/png;base64," + reag.img}
                alt=""
                style={{
                  maxWidth: "100%",
                  height: "auto",
                  display: "block",
                }}
              />
            </div>
            <div className="reagentItem">{reag.compound_class_1}</div>
          </div>
        );
      }
    }
    return <div className="selected-sequence">{arr}</div>;
  }

  function setSelectedNodeRandomly() {
    // collect node.target_found == true
    let targetNodes = network.nodes.filter((node) => node.product === true);

    setSelectedNode(
      targetNodes[Math.floor(Math.random() * targetNodes.length)]
    );
  }

  function renderCommands() {
    return (
      <div className="commands">
        <HeatmapCanvas
          data={wellToSequenceCounts}
          setter={setSelectedSequenceItems}
          sel={selectedSequenceItems}
          classMap={classMap}
        />
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
        <div className="selected">
          {renderNodeSelectedInfo()}
          {renderSelectedSequenceItems()}
        </div>
      </div>
    </div>
  );
}

// ========================================

ReactDOM.render(<Main />, document.getElementById("reactEntry"));
