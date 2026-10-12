import fs from "node:fs";

const html=fs.readFileSync("docs/index.html","utf8");
const required=["style.css","app.js","viewer","file-input","xyz-input","analyze","point-group","operations"];
for(const token of required)if(!html.includes(token))throw new Error(`Website is missing required asset or element: ${token}`);
for(const asset of ["docs/style.css","docs/app.js","docs/pgmol-core.mjs","docs/.nojekyll"]){if(!fs.existsSync(asset))throw new Error(`Website asset is missing: ${asset}`)}
if(/https?:\/\/(?!github\.com)/.test(html))throw new Error("Website unexpectedly contacts a third-party origin");
if(!html.includes("Content-Security-Policy"))throw new Error("Website has no Content Security Policy");
if(!html.includes("<noscript>"))throw new Error("Website has no JavaScript fallback message");
console.log("Static website assets and privacy constraints passed");
