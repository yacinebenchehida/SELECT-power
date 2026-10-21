// Reproduces the simulator-based numbers in SIMULATION_GUIDE.md.
// Run from the repository root:  node examples/guide_examples.js
// Uses 2,000 replicates under selection and 4,000 neutral (seed 1); takes about a minute.

const S=require('../sim.js');
const D={system:'haplodiploid',volt:1,Nf:1000,Nm:1000,noDrift:false,p0:0.01,s:0.05,h:0.5,onset:1935,years:[1935,1955,1975,2025],n:[20,20,20,40],perfect:false,depthEarly:10,depthLast:15,err:0.005,minDP:3,alpha:0.05,m:1,sMale:1,decline:1,mig:0,pMig:-1,kSib:1,maleFrac:0,reps:2000,seed:1};
const R=(o)=>{const r=S.runPower({...D,...o});return r};
const f=x=>(x*100).toFixed(1)+'%';
let r=R({}); console.log('DEFAULT hap', f(r.trend.powerEmp), f(r.two.powerEmp), 'fpr',f(r.trend.fprNaive),'het',f(r.hetAcc), 'det end', S.deterministic(D).slice(-1)[0].toFixed(3), 'Ne', S.neFor(D), 'Hlost', f(1-S.hetRetained(D)));
r=R({system:'diploid'}); console.log('DEFAULT dip', f(r.trend.powerEmp), f(r.two.powerEmp), 'det end', S.deterministic({...D,system:'diploid'}).slice(-1)[0].toFixed(3));
console.log('-- no drift vs drift (diploid, p0 .3, s .01, n 20/20/20/40, perfect)');
for(const nd of [true,false]){ r=R({system:'diploid',p0:0.3,s:0.01,perfect:true,noDrift:nd,Nf:5000,Nm:5000}); console.log(nd?'nodrift':'drift Ne10k', f(r.trend.powerEmp), f(r.two.powerEmp),'crit2',r.two.crit.toFixed(2),'fpr',f(r.trend.fprNaive)); }
console.log('-- census-like huge n with drift (diploid p0 .3 s .01 Ne 10k)');
r=R({system:'diploid',p0:0.3,s:0.01,perfect:true,Nf:5000,Nm:5000,n:[500,500,500,500]}); console.log('n500', f(r.trend.powerEmp), f(r.two.powerEmp));
r=R({system:'diploid',p0:0.3,s:0.01,perfect:true,noDrift:true,n:[500,500,500,500]}); console.log('n500 nodrift', f(r.trend.powerEmp));
console.log('-- naive FPR vs breeding numbers (hap, p0 .3, s .02)');
for(const N of [100,300,1000,3000,30000]){ r=R({p0:0.3,s:0.02,Nf:N,Nm:N}); console.log(N, 'Ne',S.neFor({system:'haplodiploid',Nf:N,Nm:N}), 'pow',f(r.trend.powerEmp),'naivePow',f(r.trend.powerNaive),'fpr',f(r.trend.fprNaive)); }
console.log('-- sample size factor (defaults hap)');
for(const k of [0.5,1,2,4]){ r=R({n:D.n.map(v=>Math.round(v*k))}); console.log(k, D.n.map(v=>Math.round(v*k)).join('/'), f(r.trend.powerEmp)); }
console.log('-- sample size factor (diploid p0 .3 s .01 Ne 10k)');
for(const k of [0.5,1,2,4]){ r=R({system:'diploid',p0:0.3,s:0.01,Nf:5000,Nm:5000,n:D.n.map(v=>Math.round(v*k))}); console.log(k, f(r.trend.powerEmp)); }
console.log('-- dominance hap vs dip');
for(const h of [0,0.5,1]) for(const sys of ['haplodiploid','diploid']){ r=R({h,system:sys}); console.log(h,sys,f(r.trend.powerEmp),'det',S.deterministic({...D,h,system:sys}).slice(-1)[0].toFixed(4)); }
console.log('-- depth');
for(const d of [2,3,6,10,20]){ r=R({depthEarly:d,depthLast:d}); console.log(d,'pow',f(r.trend.powerEmp),'het',f(r.hetAcc),'miss',f(r.missRate)); r=R({depthEarly:d,depthLast:d,maleFrac:1}); console.log('  males',f(r.trend.powerEmp),'hap',f(r.hapAcc),'miss',f(r.missRate)); }
console.log('-- ksib');
for(const k of [1,2,5,10]){ r=R({kSib:k}); console.log(k,f(r.trend.powerEmp),'fpr',f(r.trend.fprNaive)); }
console.log('-- decline'); for(const d of [1,0.3,0.1]){ r=R({decline:d}); console.log(d,f(r.trend.powerEmp),'Hlost',f(1-S.hetRetained({...D,decline:d}))); }
console.log('-- migration'); for(const m of [0,0.01,0.05]){ r=R({mig:m}); console.log(m,f(r.trend.powerEmp)); }
console.log('-- onset'); for(const y of [1935,1975,2000]){ r=R({onset:y}); console.log(y,f(r.trend.powerEmp),f(r.two.powerEmp)); }
console.log('-- m'); for(const m of [1,10,100]){ r=R({m}); console.log(m,f(r.trend.powerEmp),r.trend.method,r.trend.crit.toFixed(2)); }
console.log('-- sMale'); for(const sm of [0,1]){ r=R({sMale:sm}); console.log(sm,f(r.trend.powerEmp)); }
console.log('-- males sampling'); for(const [lab,o] of [['fem',{}],['males same ind',{maleFrac:1}],['males same alleles',{maleFrac:1,n:[40,40,40,80]}]]){ r=R(o); console.log(lab,f(r.trend.powerEmp)); }
console.log('-- 4pt vs 2pt late onset'); r=R({onset:1980}); console.log(f(r.trend.powerEmp),f(r.two.powerEmp));
