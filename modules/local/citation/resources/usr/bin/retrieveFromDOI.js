#!/usr/bin/env node

// Load Citation.js
const { Cite } = require("@citation-js/core");
// Load plugins
require("@citation-js/plugin-doi");
require("@citation-js/plugin-csl");
const TurndownService = require("turndown");

// Parse CLI arguments
const args = process.argv.slice(2);
const usage = `
Usage: node main.js <DOI> [--format <format>] [--template <template>]

Options:
  --format    Output format (html, markdown; default: html)
  --template  Citation template (default: apa)
`;

if (args.length === 0 || args.includes("--help")) {
  console.log(usage);
  process.exit(0);
}

const doi = args[0];
let format = "html";
let template = "apa";

for (let i = 1; i < args.length; i++) {
  if (args[i] === "--format" && args[i + 1]) {
    format = args[i + 1];
    i++;
  } else if (args[i] === "--template" && args[i + 1]) {
    template = args[i + 1];
    i++;
  }
}

async function generateCitation(doi, format = "html", template = "apa") {
  try {
    const exampleCite = await Cite.async(doi);
    let formatCite = html;
    if (format === "markdown") {
      formatCite = "html";
    } else {
      formatCite = format;
    }
    const output = exampleCite.format("bibliography", {
      format: formatCite,
      template,
    });
    if (format === "markdown") {
      const turndownService = new TurndownService();
      console.log(turndownService.turndown(output));
    } else {
      console.log(output);
    }
  } catch (error) {
    console.error("Failed to retrieve citation data:", error);
  }
}

// Example with a specific DOI
generateCitation(doi, format, template);
