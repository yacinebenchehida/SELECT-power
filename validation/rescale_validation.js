// Validation of SPOT's rescaling: full (exact) simulations vs rescaled ones on the same designs.
// Run from the repository root: node validation/rescale_validation.js > validation/rescale_results.tsv
const S=require('../sim.js');
const base={system:'haploid',genTime:1,Nf:1e4,Nm:1e4,noDrift:false,p0:0.1,s:0.0005,h:0.5,onset:0,n:[20,20,20,40],perfect:false,depthEarly:10,depthLast:15,err:0.005,minDP:3,alpha:0.05,m:1,sMale:1,decline:1,mig:0,pMig:-1,kSib:1,maleFrac:0,reps:1000,seed:1};
const designFor=G=>({years:[0,Math.round(G/3),Math.round(2*G/3),G]});
// distribution of the true population frequency at the last sampling point (no sampling), exact vs rescaled
function finalFreq(P,s,reps,seed){ const Q=S.prep(Object.assign({},P,{n:P.n.map(()=>0)})), R=S.makeRng(seed); const out=new Float64Array(reps);
  for(let r=0;r<reps;r++){ const o=S.simulateRep(Q,s,R,true); out[r]=o.traj[o.traj.length-1]; } return out; }
const stats=a=>{let m=0;for(const x of a)m+=x;m/=a.length;let v=0;for(const x of a)v+=(x-m)*(x-m);v/=a.length-1;let lost=0,fix=0;for(const x of a){if(x<=0)lost++;if(x>=1)fix++}return {mean:m,sd:Math.sqrt(v),lost:lost/a.length,fixed:fix/a.length}};
const rows=[['G','system','N','s','p0','lambda','scaled_s','scaled_N','power_exact','power_rescaled','two_exact','two_rescaled','crit_exact','crit_rescaled','naiveFPR_exact','naiveFPR_rescaled','meanFinal_exact','meanFinal_rescaled','sdFinal_exact','sdFinal_rescaled','neutralSD_exact','neutralSD_rescaled','neutralSD_theory','ms_exact','ms_rescaled']];
const t0=Date.now();
for(const G of [5000,20000,60000]) for(const system of ['haploid','diploid']) for(const N of [1e3,1e4,1e5]) for(const s of [0.0001,0.0005,0.002]) for(const p0 of [0.01,0.1]){
  if(G===60000&&(N<1e4||s>0.0005||p0<0.1)) continue; // long runs: a subset, to keep the exact simulations affordable
  const P=Object.assign({},base,designFor(G),{system,Nf:N,Nm:N,s,p0});
  const sc=S.rescale(P); if(sc.lambda===1) continue;
  let t=Date.now(); const e=S.runPower(Object.assign({},P,{noRescale:true})); const te=Date.now()-t;
  t=Date.now(); const r=S.runPower(Object.assign({},P,{seed:2})); const tr=Date.now()-t;
  const fe=stats(finalFreq(Object.assign({},P,{noRescale:true}),s,1000,3)), fr=stats(finalFreq(sc.P,sc.P.s,1000,4));
  const ne=stats(finalFreq(Object.assign({},P,{noRescale:true}),0,1000,5)), nr=stats(finalFreq(sc.P,0,1000,6));
  const Ne=S.neFor(P), loss=system==='haploid'?1-Math.pow(1-1/Ne,G):1-Math.pow(1-1/(2*Ne),G), sdTheory=Math.sqrt(p0*(1-p0)*loss);
  rows.push([G,system,N,s,p0,sc.lambda,+sc.P.s.toPrecision(3),+sc.P.Nf.toPrecision(3),e.trend.powerEmp,r.trend.powerEmp,e.two.powerEmp,r.two.powerEmp,+e.trend.crit.toFixed(3),+r.trend.crit.toFixed(3),e.trend.fprNaive,r.trend.fprNaive,
    +fe.mean.toFixed(4),+fr.mean.toFixed(4),+fe.sd.toFixed(4),+fr.sd.toFixed(4),+ne.sd.toFixed(4),+nr.sd.toFixed(4),+sdTheory.toFixed(4),te,tr]);
  console.error(rows.length-1, G,system,N,s,p0,'lambda',sc.lambda,'power',e.trend.powerEmp,r.trend.powerEmp,((Date.now()-t0)/1000).toFixed(0)+'s');
}
console.log(rows.map(r=>r.join('\t')).join('\n'));
