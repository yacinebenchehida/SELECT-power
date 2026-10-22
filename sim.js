/* ===== SPOT (Selection Power Over Time): core (pure functions, no DOM) =====
   Usable in the browser or in Node (module.exports at the bottom). See README.md. */
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296}}
function makeRng(seed){
  const u=mulberry32(seed>>>0); let spare=null;
  const gauss=()=>{if(spare!==null){const v=spare;spare=null;return v}let a,b,r;do{a=2*u()-1;b=2*u()-1;r=a*a+b*b}while(r>=1||r===0);const f=Math.sqrt(-2*Math.log(r)/r);spare=b*f;return a*f};
  function binom(n,p){
    if(n<=0||p<=0) return 0; if(p>=1) return n;
    if(p>0.5) return n-binom(n,1-p);
    if(n<=40){let k=0;for(let i=0;i<n;i++) if(u()<p) k++;return k}
    const m=n*p;
    if(m<20){ // inversion
      const q=1-p, r=p/q; let f=Math.pow(q,n), c=f, k=0, x=u();
      while(x>c&&k<n){f*=r*(n-k)/(k+1);k++;c+=f}
      return k;
    }
    const k=Math.round(m+Math.sqrt(m*(1-p))*gauss());
    return k<0?0:k>n?n:k;
  }
  function pois(l){
    if(l<=0) return 0;
    if(l>60){const k=Math.round(l+Math.sqrt(l)*gauss());return k<0?0:k}
    const L=Math.exp(-l);let k=0,p=1;do{k++;p*=u()}while(p>L);return k-1;
  }
  return {u,gauss,binom,pois};
}
const YEARS=[1935,1955,1975,2025]; // default sampling years; override with P.years (any number of points, >= 2)
/* years per generation: P.genTime, or 1/P.volt (generations per year) for older parameter sets */
function genTimeOf(P){return P.genTime>0?P.genTime:1/(P.volt||1)}
/* generation index of each sampling year */
function gensOf(P){const Yr=P.years||YEARS, gt=genTimeOf(P); return Yr.map(y=>Math.round((y-Yr[0])/gt))}

/* breeding-number multiplier in generation g: per-interval Ne (P.neByBlock) or exponential change (P.decline) */
function sizeFactor(P,gens,g){
  const K=gens.length, G=gens[K-1];
  if(P.neByBlock&&P.neBlocks&&P.neBlocks.length>=K-1){let j=0; while(j<K-2&&g>=gens[j+1]) j++; return P.neBlocks[j]/neFor(P)}
  return Math.pow(P.decline,g/G);
}
/* one replicate. Returns allele counts at the sampled points + optional trajectories (female, male) */
function simulateRep(P,s,R,keep){
  const Yr=P.years||YEARS, gens=gensOf(P), K=gens.length, G=gens[K-1];
  const onset=(P.onset-Yr[0])/genTimeOf(P), hap=P.system==='haplodiploid', hpl=P.system==='haploid', h=P.h, sMrel=P.sMale;
  const inf=P.N===Infinity, mig=P.mig, pMig=P.pMig;
  let pf=P.p0, pm=P.p0, si=0;
  const x=new Array(K).fill(0), n=new Array(K).fill(0), traj=keep?[]:null, trajM=keep?[]:null;
  let hetTrue=0, hetCalled=0, missing=0, sampled=0, hapTot=0, hapOk=0;
  const mf=P.maleFrac||0;
  const callOne=(gt,dep)=>{ // returns called genotype or -1 (missing)
    sampled++;
    if(P.perfect){ if(gt===1){hetTrue++;hetCalled++} return gt; }
    const d=R.pois(dep);
    if(d<P.minDP){missing++;return -1}
    if(gt===1) hetTrue++;
    const q=gt===0?P.err:gt===1?0.5:1-P.err;
    const k=R.binom(d,q);
    const l0=k*P.lE+(d-k)*P.l1E, l1=-d*Math.LN2, l2=k*P.l1E+(d-k)*P.lE;
    const c=(l1>=l0&&l1>=l2)?1:(l2>l0?2:0);
    if(gt===1&&c===1) hetCalled++;
    return c;
  };
  const callHap=(a,dep)=>{ // haploid male: returns called allele or -1
    sampled++;
    if(P.perfect){hapTot++;hapOk++;return a}
    const d=R.pois(dep);
    if(d<P.minDP){missing++;return -1}
    hapTot++;
    const k=R.binom(d,a?1-P.err:P.err);
    const c=2*k>d?1:(2*k<d?0:(R.u()<0.5?1:0));
    if(c===a) hapOk++;
    return c;
  };
  for(let g=0; g<=G; g++){
    const sg=g>onset?s:0, sm=sg*sMrel;
    let pfs, pms, fAA=0, fAa=0, mAA=0, mAa=0;
    if(hpl){ pfs=pf*(1+sg)/(1+pf*sg); pms=pfs; } // haploid: one allele per individual, fitness 1+s or 1
    else{
      const AA=pf*pm, Aa=pf*(1-pm)+pm*(1-pf);
      const w=AA*(1+sg)+Aa*(1+h*sg)+(1-AA-Aa);
      fAA=AA*(1+sg)/w; fAa=Aa*(1+h*sg)/w; pfs=fAA+fAa/2;
      if(hap){ pms=pf*(1+sm)/(pf*(1+sm)+1-pf); }
      else { const wm=AA*(1+sm)+Aa*(1+h*sm)+(1-AA-Aa); mAA=AA*(1+sm)/wm; mAa=Aa*(1+h*sm)/wm; pms=mAA+mAa/2; }
    }
    if(keep){traj.push(pfs);trajM.push(pms)}
    while(si<K&&g===gens[si]){ // several sampling years can fall in the same generation
      const ni=P.n[si], dep=si===K-1?P.depthLast:P.depthEarly, k=P.kSib;
      if(hpl){
        for(let i=0;i<ni;i++){const c=callHap(R.u()<pfs?1:0,dep);if(c>=0){x[si]+=c;n[si]+=1}}
      }else if(k<=1){
        for(let i=0;i<ni;i++){
          if(mf>0&&R.u()<mf){ // male
            if(hap){const c=callHap(R.u()<pms?1:0,dep);if(c>=0){x[si]+=c;n[si]+=1}}
            else{const r=R.u();const gt=r<mAA?2:(r<mAA+mAa?1:0);const c=callOne(gt,dep);if(c>=0){x[si]+=c;n[si]+=2}}
          }else{const r=R.u();const gt=r<fAA?2:(r<fAA+fAa?1:0);const c=callOne(gt,dep);if(c>=0){x[si]+=c;n[si]+=2}}
        }
      }else{
        let left=ni;
        while(left>0){
          const r=R.u(); const mg=r<fAA?2:(r<fAA+fAa?1:0);         // mother (breeding female)
          let fa, fg=0;                                              // father
          if(hap) fa=R.u()<pms?1:0; else {const q=R.u(); fg=q<mAA?2:(q<mAA+mAa?1:0);}
          const kk=Math.min(k,left); left-=kk;
          for(let j=0;j<kk;j++){
            const mat=mg===2?1:mg===0?0:(R.u()<0.5?1:0);
            const pat=hap?fa:(fg===2?1:fg===0?0:(R.u()<0.5?1:0));
            if(hap&&mf>0&&R.u()<mf){const c=callHap(mat,dep);if(c>=0){x[si]+=c;n[si]+=1};continue} // son: maternal allele only
            const c=callOne(mat+pat,dep); if(c>=0){x[si]+=c;n[si]+=2}
          }
        }
      }
      si++;
    }
    if(inf){pf=pfs;pm=pms}
    else{
      const f=sizeFactor(P,gens,g), Nf=Math.max(2,Math.round(P.Nf*f)), Nm=Math.max(2,Math.round(P.Nm*f));
      if(hpl){ pf=R.binom(Nf,pfs)/Nf; pm=pf; }
      else{
        pf=R.binom(2*Nf,pfs)/(2*Nf);
        pm=hap? R.binom(Nm,pms)/Nm : R.binom(2*Nm,pms)/(2*Nm);
      }
    }
    if(mig>0){const pmg=pMig<0?P.p0:pMig; pf=(1-mig)*pf+mig*pmg; pm=(1-mig)*pm+mig*pmg;}
  }
  return {x,n,traj,trajM,hetTrue,hetCalled,missing,sampled,hapTot,hapOk};
}
/* Cochran–Armitage trend Z on allele counts; idx = which points to use */
function trendZ(x,n,t,idx){
  let N=0,X=0,T=0;
  for(const i of idx){N+=n[i];X+=x[i];T+=n[i]*t[i]}
  if(N===0) return 0;
  const pb=X/N, tb=T/N; if(pb<=0||pb>=1) return 0;
  let num=0,den=0;
  for(const i of idx){num+=x[i]*(t[i]-tb);den+=n[i]*(t[i]-tb)*(t[i]-tb)}
  if(den<=0) return 0;
  return num/Math.sqrt(pb*(1-pb)*den);
}
/* inverse normal CDF (Acklam) */
function qnorm(p){
  const a=[-39.69683028665376,220.9460984245205,-275.9285104469687,138.357751867269,-30.66479806614716,2.506628277459239],
  b=[-54.47609879822406,161.5858368580409,-155.6989798598866,66.80131188771972,-13.28068155288572],
  c=[-.007784894002430293,-.3223964580411365,-2.400758277161838,-2.549732539343734,4.374664141464968,2.938163982698783],
  d=[.007784695709041462,.3224671290700398,2.445134137142996,3.754408661907416];
  const pl=.02425; let q,r;
  if(p<pl){q=Math.sqrt(-2*Math.log(p));return(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)}
  if(p<=1-pl){q=p-.5;r=q*q;return(((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)}
  q=Math.sqrt(-2*Math.log(1-p));return-(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
}
/* ---------- rescaling for very many generations ----------
   When the sampling period spans more than SIM_TARGET generations, every simulated generation stands for
   lambda real ones: selection becomes (1+s)^lambda - 1, breeding numbers (and Ne per interval) are divided by
   lambda, migration becomes 1-(1-m)^lambda and sampling times are divided by lambda. This keeps Ne*s and the
   total drift over the period unchanged (the standard rescaling used in forward simulators such as SLiM).
   lambda is capped so that scaled selection stays <= MAX_SCALED_S and scaled breeding numbers >= MIN_SCALED_N
   (limits set from validation/rescale_limits.js);
   if those caps would leave more than SIM_HARD generations, lambda is raised anyway and flagged. */
const SIM_TARGET=2000, SIM_HARD=10000, MAX_SCALED_S=0.025, MIN_SCALED_N=200;
function rescale(P){
  const g=gensOf(P), G=g[g.length-1];
  if(P.noRescale||G<=SIM_TARGET) return {P,lambda:1,G,Gsim:G,outside:false};
  const want=Math.ceil(G/SIM_TARGET), need=Math.ceil(G/SIM_HARD);
  let lam=want;
  const sAbs=Math.max(Math.abs(P.s||0),Math.abs((P.s||0)*(P.sMale==null?1:P.sMale)));
  if(sAbs>0) lam=Math.min(lam,Math.max(1,Math.floor(Math.log(1+MAX_SCALED_S)/Math.log(1+sAbs))));
  if(!P.noDrift){
    const nMin=(P.system==='haploid'?P.Nf:Math.min(P.Nf,P.Nm))*(P.neByBlock&&P.neBlocks?Math.min(...P.neBlocks)/neFor(P):Math.min(1,P.decline));
    lam=Math.min(lam,Math.max(1,Math.floor(nMin/MIN_SCALED_N)));
  }
  if(P.mig>0&&P.mig<1) lam=Math.min(lam,Math.max(1,Math.floor(Math.log(0.5)/Math.log(1-P.mig))));
  let outside=false; if(lam<need){lam=need; outside=true;}
  if(lam<=1) return {P,lambda:1,G,Gsim:G,outside:false};
  const Q=Object.assign({},P); delete Q.volt;
  Q.genTime=genTimeOf(P)*lam;
  Q.s=Math.pow(1+P.s,lam)-1;
  const sm=P.s*(P.sMale==null?1:P.sMale); Q.sMale=Q.s>0?(Math.pow(1+sm,lam)-1)/Q.s:P.sMale;
  Q.Nf=P.Nf/lam; Q.Nm=P.Nm/lam; if(P.neBlocks) Q.neBlocks=P.neBlocks.map(v=>v/lam);
  Q.mig=1-Math.pow(1-P.mig,lam);
  const gs=gensOf(Q);
  return {P:Q,lambda:lam,G,Gsim:gs[gs.length-1],outside};
}
function prep(P){
  const Q=Object.assign({},P);
  Q.lE=Math.log(Q.err); Q.l1E=Math.log(1-Q.err);
  if(Q.noDrift) Q.N=Infinity; else Q.N=Q.Nf;
  return Q;
}
/* run R selected + 2R neutral replicates; compute power several ways */
function runPower(P0,opt={}){
  const sc=rescale(P0), P=prep(sc.P), R=makeRng(P.seed), reps=P.reps, nn=reps*2;
  const t=gensOf(P), K=t.length;
  const ALL=t.map((_,i)=>i), TWO=[0,K-1];
  const keepK=opt.keep||0;
  const sel={trend:new Float64Array(reps),two:new Float64Array(reps)}, neu={trend:new Float64Array(nn),two:new Float64Array(nn)};
  const trajS=[],trajN=[]; let het=0,hetOk=0,miss=0,samp=0,hT=0,hO=0;
  for(let r=0;r<reps;r++){const o=simulateRep(P,P.s,R,r<keepK);sel.trend[r]=trendZ(o.x,o.n,t,ALL);sel.two[r]=trendZ(o.x,o.n,t,TWO);if(o.traj)trajS.push(o.traj);het+=o.hetTrue;hetOk+=o.hetCalled;miss+=o.missing;samp+=o.sampled;hT+=o.hapTot;hO+=o.hapOk}
  for(let r=0;r<nn;r++){const o=simulateRep(P,0,R,r<keepK);neu.trend[r]=trendZ(o.x,o.n,t,ALL);neu.two[r]=trendZ(o.x,o.n,t,TWO);if(o.traj)trajN.push(o.traj)}
  const a=P.alpha/P.m, zc=qnorm(1-a);
  const res={alphaAdj:a,zNaive:zc,trajS,trajN,sel,neu,lambda:sc.lambda,genReal:sc.G,genSim:sc.Gsim,rescaleOutside:sc.outside,
    hetAcc:P.perfect?1:(het?hetOk/het:NaN), missRate:P.perfect?0:miss/samp, hapAcc:hT?hO/hT:NaN};
  for(const k of ['trend','two']){
    const ns=Float64Array.from(neu[k]).sort();
    let crit, method;
    if(a*nn>=10){crit=ns[Math.min(nn-1,Math.ceil((1-a)*nn)-1)];method='empirical'}
    else{let m=0,v=0;for(const z of ns)m+=z;m/=nn;for(const z of ns)v+=(z-m)*(z-m);v/=nn-1;crit=m+Math.sqrt(v)*zc;method='gaussian'}
    let pe=0,pn=0,fp=0;
    for(const z of sel[k]){if(z>crit)pe++;if(z>zc)pn++}
    for(const z of neu[k]) if(z>zc) fp++;
    res[k]={crit,method,powerEmp:pe/reps,powerNaive:pn/reps,fprNaive:fp/nn};
  }
  return res;
}
function deterministic(P0,both){ // exact up to 200,000 generations, rescaled beyond
  const g=gensOf(P0), D=Object.assign({},P0,{noDrift:true,perfect:true});
  const P=prep(g[g.length-1]>200000?rescale(D).P:D);
  const R=makeRng(1); const o=simulateRep(Object.assign(P,{n:(P.years||YEARS).map(()=>0)}),P.s,R,true); return both?{f:o.traj,m:o.trajM}:o.traj;
}
/* neutral diversity retained between first and last sampling point: prod (1 - 1/(2 Ne_g)), or prod (1 - 1/N_g) for haploids */
function hetRetained(P){
  if(P.noDrift) return 1;
  const gens=gensOf(P), G=gens[gens.length-1], step=Math.max(1,Math.ceil(G/100000)); let H=1;
  for(let g=0;g<G;g+=step){const k=Math.min(step,G-g), f=sizeFactor(P,gens,g), ne=neFor({system:P.system,Nf:Math.max(2,P.Nf*f),Nm:Math.max(2,P.Nm*f)}); H*=Math.pow(P.system==='haploid'?1-1/ne:1-1/(2*ne),k)}
  return H;
}
/* Ne: diploid 4NfNm/(Nf+Nm); haplodiploid 9NfNm/(4Nm+2Nf); haploid N (= Nf, number of breeding individuals) */
function neFor(P){return P.system==='haploid'? P.Nf : P.system==='haplodiploid'? 9*P.Nf*P.Nm/(4*P.Nm+2*P.Nf) : 4*P.Nf*P.Nm/(P.Nf+P.Nm)}
if(typeof module!=='undefined') module.exports={runPower,deterministic,neFor,hetRetained,simulateRep,makeRng,prep,qnorm,gensOf,genTimeOf,sizeFactor,rescale,YEARS};
