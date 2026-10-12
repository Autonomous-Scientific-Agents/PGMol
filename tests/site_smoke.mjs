import fs from "node:fs";
import path from "node:path";

const html=fs.readFileSync("docs/index.html","utf8");
const required=["style.css","app.js","viewer","file-input","xyz-input","analyze","point-group","operations"];
for(const token of required)if(!html.includes(token))throw new Error(`Website is missing required asset or element: ${token}`);
for(const asset of ["docs/style.css","docs/app.js","docs/pgmol-core.mjs","docs/.nojekyll","docs/agent-traces/index.html","docs/agent-traces/pgmol-point-group-analyzer/replay.html"]){if(!fs.existsSync(asset))throw new Error(`Website asset is missing: ${asset}`)}
if(!html.includes("agent-traces/"))throw new Error("Website does not link to the agent traces");
if(/https?:\/\/(?!github\.com)/.test(html))throw new Error("Website unexpectedly contacts a third-party origin");
if(!html.includes("Content-Security-Policy"))throw new Error("Website has no Content Security Policy");
if(!html.includes("<noscript>"))throw new Error("Website has no JavaScript fallback message");

function walk(directory){return fs.readdirSync(directory,{withFileTypes:true}).flatMap(entry=>{const file=path.join(directory,entry.name);return entry.isDirectory()?walk(file):[file]})}
const published=walk("docs/agent-traces");
if(published.some(file=>file.endsWith(".json")||file.endsWith(".jsonl")))throw new Error("Raw trace source must not be published through Pages");
for(const file of published.filter(file=>file.endsWith(".html"))){
  const page=fs.readFileSync(file,"utf8");
  if(!page.includes("Content-Security-Policy"))throw new Error(`${file} has no Content Security Policy`);
  const external=[...page.matchAll(/https?:\/\/[^\s"'<>]+/g)].map(match=>match[0]).filter(url=>url!=="http://www.w3.org/2000/svg"&&!(file.endsWith("index.html")&&url.startsWith("https://github.com/")));
  if(external.length)throw new Error(`${file} unexpectedly references external origins: ${external.join(", ")}`);
}
console.log("Static website assets, trace publication, and privacy constraints passed");
