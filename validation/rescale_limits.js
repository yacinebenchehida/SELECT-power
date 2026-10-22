// Stress test: force increasing scaling factors and compare with the exact simulation.
// Run from the repository root: node validation/rescale_limits.js
const S=require('../sim.js');
const base={system:'haploid',genTime:1,Nf:1e4,Nm:1e4,noDrift:false,p0:0.1,s:0.0005,h:0.5,onset:0,years:[0,6667,13333,20000],n:[20,20,20,40],perfect:false,depthEarly:10,depthLast:15,err:0.005,minDP:3,alpha:0.05,m:1,sMale:1,decline:1,mig:0,pMig:-1,kSib:1,maleFrac:0,reps:2000,seed:1};
function force(P,lam){const Q=Object.assign({},P,{noRescale:true}); Q.genTime=P.genTime*lam; Q.s=Math.pow(1+P.s,lam)-1; Q.Nf=P.Nf/lam; Q.Nm=P.Nm/lam; return Q;}
const cases=[['haploid N=1e4 s=5e-4 p0=0.1',{}],['haploid N=1e3 s=2e-3 p0=0.1',{Nf:1e3,Nm:1e3,s:0.002}],['diploid N=1e5 s=5e-4 p0=0.01',{system:'diploid',Nf:1e5,Nm:1e5,p0:0.01}],['diploid N=1e3 s=5e-4 p0=0.1',{system:'diploid',Nf:1e3,Nm:1e3}]];
for(const [lab,o] of cases){ const P=Object.assign({},base,o); const ex=S.runPower(Object.assign({},P,{noRescale:true})).trend.powerEmp; const out=[];
  for(const lam of [2,5,10,20,50,100,200]){ const Q=force(P,lam); if(Q.Nf<2) break; const r=S.runPower(Object.assign({},Q,{seed:2})).trend.powerEmp; out.push(`λ=${lam} (s'=${Q.s.toPrecision(2)}, N'=${Math.round(Q.Nf)}): ${(r*100).toFixed(1)}`); }
  console.log(`${lab} | exact ${(ex*100).toFixed(1)} | `+out.join(' | ')); }
