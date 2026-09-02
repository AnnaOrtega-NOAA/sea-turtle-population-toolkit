// Expected-value teaching kernels, transcribed from the inspected Methods script (8)
// and app 2026-09-02. No model fitting is performed in the browser.
export const params = Object.freeze({cf:5.5,eggs:77.9,ri:3.06,linf:142.7,k:.2262,t0:-.17,lmat:139.1325,sig_mat:6.3399,pf:.73,pj:.81,pa:.893,mat_p:.99,year1:.0625});
export const sum = a => a.reduce((s,v)=>s+v,0);
export const cumsum = a => {let s=0;return a.map(v=>s+=v)};
export const transition = (l,p=params) => l>=p.mat_p*p.linf ? 1 : 1/(1+Math.exp(-(l-p.lmat)/p.sig_mat));
export const ageFromLength=(l,p=params)=>p.t0-Math.log(1-l/p.linf)/p.k;
export const lengthAt=(age,p=params)=>p.linf*(1-Math.exp(-p.k*(age-p.t0)));
export function takePath(length,p=params,horizon=60,mortality=1,count=1){
 let age0=ageFromLength(length,p), maxAge=ageFromLength(.99*p.linf,p); if(!Number.isFinite(age0))age0=maxAge;
 const atMax=age0>=maxAge || Math.abs(age0-maxAge)<1e-12;
 const growthCount=Math.floor(maxAge-age0+1e-12)+1;
 let survival=1, notNested=1;
 const rows=Array.from({length:horizon},(_,j)=>{
 const l=(atMax || (growthCount<horizon&&j>=growthCount))?p.linf*p.mat_p:lengthAt(age0+j,p);
 const pt=transition(l,p); survival*=((1-pt)*p.pj+pt*p.pa);const fn=pt*notNested;notNested*=1-pt;
 return {year:j+1,length:l,pt,fn,survival,ane:count*survival*fn*p.pf/p.ri*mortality};});
 return rows;
}
export function givePath(type,amount,p=params,horizon=60,releaseAge=1){
 const maxAge=ageFromLength(p.mat_p*p.linf,p),age0=type==='adult'?maxAge:type==='nest'?0:releaseAge;
 const len=age=>age>=maxAge?p.linf*p.mat_p:lengthAt(age,p);
 let survival=1, notNested=1;
 return Array.from({length:horizon},(_,j)=>{const intervalPt=transition(len(age0+j),p);let phi=(1-intervalPt)*p.pj+intervalPt*p.pa;if(j===0&&type==='nest')phi=p.year1;
 survival*=phi;const pt=transition(len(age0+j+1),p),fn=pt*notNested;notNested*=1-pt;
 return {year:j+1,fn,survival,ane:amount*survival*fn*p.pf/p.ri};});
}
export function summary(rows){const a=rows.map(d=>d.ane),total=sum(a);let index=cumsum(a).findIndex(d=>d>=total*.5);return {total,survive:sum(rows.map(d=>d.survival*d.fn)),median:total>0?index+1:null};}
export const step=(n,t,g,u,q,z)=>Math.max(0,Math.max(0,n-t)*Math.exp(u)+Math.sqrt(q)*z+g);
export function rng(seed){return ()=>{seed=(Math.imul(1664525,seed)+1013904223)>>>0;return(seed+.5)/4294967296}}
export function future(seed=20260827,take=3,give=2,n=160,horizon=30){
 const random=rng(seed),normal=()=>Math.sqrt(-2*Math.log(random()))*Math.cos(2*Math.PI*random());
 // Synthetic posterior-like rows solely for illustrating pairing; never fitted results.
 const posterior=Array.from({length:256},(_,i)=>{let v=normal();return {id:i+1,N:300*Math.exp(.12*normal()),U:-.015+.016*v,Q:.0025*Math.exp(.45*v)}});
 const all=[[],[],[],[]];let ledger;
 for(let s=0;s<n;s++){const start=Math.floor(random()*posterior.length),rest=posterior.map((_,i)=>i);for(let j=rest.length-1;j>0;j--){const k=Math.floor(random()*(j+1));[rest[k],rest[j]]=[rest[j],rest[k]]}
 const rows=[start,...rest.slice(0,horizon-1)].map((i,j)=>({...posterior[i],year:j+1,z:normal()}));
 if(s===0)ledger=rows;
 const paths=Array.from({length:4},()=>[posterior[start].N]);
 rows.forEach(r=>paths.forEach((path,k)=>path.push(step(path.at(-1),k===1||k===3?take:0,k===2||k===3?give:0,r.U,r.Q,r.z))));
 paths.forEach((p,k)=>all[k].push(p));}
 const quant=(values,p)=>{const a=[...values].sort((x,y)=>x-y);const idx=(a.length-1)*p,lo=Math.floor(idx);return a[lo]+(a[Math.ceil(idx)]-a[lo])*(idx-lo)};
 const bands=all.map(paths=>Array.from({length:horizon+1},(_,j)=>({lo:quant(paths.map(p=>p[j]),.05),mid:quant(paths.map(p=>p[j]),.5),hi:quant(paths.map(p=>p[j]),.95)})));
 const diffs=all.map(paths=>Array.from({length:horizon+1},(_,j)=>{const v=paths.map((p,s)=>p[j]-all[0][s][j]);return {lo:quant(v,.05),mid:quant(v,.5),hi:quant(v,.95)}}));
 return {all,bands,diffs,ledger};
}
