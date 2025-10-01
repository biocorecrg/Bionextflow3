#!/usr/bin/env node

// Load Citation.js
const { Cite } = require("@citation-js/core");
// Load plugins
require("@citation-js/plugin-doi");
require("@citation-js/plugin-csl");
const TurndownService = require("turndown");
const fs = require("fs");
const yaml = require("js-yaml");

// Parse CLI arguments
const args = process.argv.slice(2);
const usage = `
Usage: node addCitationFromYaml.js [--input <yaml>] [--output <yaml>] [--template <template>]

Options:
  --input     Yaml input file
  --output    Yaml output file
  --template  Citation template (default: apa)
  --format    Format for the citation (default: markdown)
`;

if (args.length === 0 || args.includes("--help")) {
  console.log(usage);
  process.exit(0);
}

// Helper to get argument value
function getArgValue(flag) {
  const idx = args.indexOf(flag);
  return idx !== -1 && args[idx + 1] ? args[idx + 1] : null;
}

const inputFile = getArgValue("--input");
const outputFile = getArgValue("--output");
let template = getArgValue("--template") || "apa";
let format = getArgValue("--format") || "markdown";

if (!inputFile || !outputFile) {
  console.error("Error: --input and --output arguments are required.");
  console.log(usage);
  process.exit(1);
}

async function addCitationsToYaml(
  inputPath,
  outputPath,
  template = "apa",
  format = "markdown",
) {
  try {
    // Read and parse YAML
    const fileContent = fs.readFileSync(inputPath, "utf8");
    const data = yaml.load(fileContent);

    async function processSection(section) {
      if (!section || typeof section !== "object") return;
      for (const key of Object.keys(section)) {
        const entry = section[key];
        if (entry && entry.doi && entry.doi.trim()) {
          try {
            const citeObj = await Cite.async(entry.doi);
            let citation;
            if (format === "markdown") {
              const htmlCitation = citeObj.format("bibliography", {
                format: "html",
                template,
              });
              const turndownService = new TurndownService();
              citation = turndownService.turndown(htmlCitation);
            } else {
              citation = citeObj.format("bibliography", {
                format,
                template,
              });
            }
            entry.citation = citation;
          } catch (e) {
            entry.citation = `Error generating citation: ${e.message}`;
          }
        }
      }
    }

    // Process all top-level sections
    for (const sectionName of Object.keys(data)) {
      await processSection(data[sectionName]);
    }

    // Write updated YAML
    const newYaml = yaml.dump(data, { lineWidth: 120 });
    fs.writeFileSync(outputPath, newYaml, "utf8");
    console.log(`Citations added and written to ${outputPath}`);
  } catch (err) {
    console.error("Failed to process YAML:", err);
    process.exit(1);
  }
}

addCitationsToYaml(inputFile, outputFile, template, format);
