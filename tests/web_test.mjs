import fs from "node:fs";
import path from "node:path";
import {analyzeMolecule,parseXYZ} from "../docs/pgmol-core.mjs";

const roots=["core_xyz","challenge_xyz"].filter(root=>fs.existsSync(root));
let passed=0;
for(const root of roots){
  for(const name of fs.readdirSync(root).filter(name=>name.endsWith(".xyz")).sort()){
    const expected=path.basename(name,".xyz").split("_").at(-1);
    const parsed=parseXYZ(fs.readFileSync(path.join(root,name),"utf8"));
    const actual=analyzeMolecule(parsed.symbols,parsed.coordinates,.1).group;
    if(actual!==expected)throw new Error(`${root}/${name}: expected ${expected}, got ${actual}`);
    passed++;
  }
}
const malformed=["", "x\ncomment\n", "1\ncomment\nXx 0 0 0\n", "2\ncomment\nH 0 0 0\n"];
for(const text of malformed){let rejected=false;try{parseXYZ(text)}catch{rejected=true}if(!rejected)throw new Error(`Malformed XYZ was accepted: ${JSON.stringify(text)}`)}
const water=parseXYZ(fs.readFileSync("core_xyz/H2O_C2v.xyz","utf8"));
const moved=water.coordinates.map(([x,y,z])=>[-y+8,x-4,z+2]);
if(analyzeMolecule([...water.symbols].reverse(),[...moved].reverse(),.1).group!=="C2v")throw new Error("Browser classifier is not rigid-motion/order invariant");
const atom=analyzeMolecule(["Ar"],[[3,-2,1]],.1);
if(atom.group!=="Kh"||atom.operationCount!=="∞")throw new Error("Single atom metadata is incorrect");
const h2=analyzeMolecule(["H","H"],[[0,0,-.37],[0,0,.37]],.1);
if(h2.group!=="Dinfh"||h2.operationCount!=="∞")throw new Error("Linear group metadata is incorrect");
const d3d=parseXYZ(fs.readFileSync("core_xyz/C2H6_D3d.xyz","utf8"));
const d3dResult=analyzeMolecule(d3d.symbols,d3d.coordinates,.1);
if(!d3dResult.operations.some(op=>op.symbol==="S6"&&op.detail.includes("60°")))throw new Error("S6 operation metadata is incorrect");
const ih=parseXYZ(fs.readFileSync("core_xyz/B12H12_Ih.xyz","utf8"));
const ihResult=analyzeMolecule(ih.symbols,ih.coordinates,.1);
if(!ihResult.operations.some(op=>op.symbol==="C5"&&op.detail.includes("144°")))throw new Error("Higher-power rotation metadata is incorrect");
console.log(`${passed} browser classifier fixtures passed; ${malformed.length} malformed inputs rejected`);
