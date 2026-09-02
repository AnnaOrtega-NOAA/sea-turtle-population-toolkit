import './site.js';

// Optional guided tour. All teaching changes are restored when the tour closes.
const $ = s => document.querySelector(s);
const all = s => [...document.querySelectorAll(s)];
if ($('#map') && !$('#start-science-tour')) {
  const css = document.createElement('link');
  css.rel = 'stylesheet'; css.href = new URL('./tour.css', import.meta.url).href;
  document.head.append(css);
  const start = document.createElement('button');
  start.id = 'start-science-tour'; start.className = 'tour-start';
  start.textContent = 'Start guided tour'; start.type = 'button';
  $('.intro').insertAdjacentElement('afterend', start);
  const click = s => $(s)?.click();
  const set = (id, value) => { const e = $('#'+id); if(e) {e.value=value;e.dispatchEvent(new Event('input',{bubbles:true}));} };
  const step = (panel, target, title, text, action) => ({panel,target,title,text,action});
  const steps = [
    step('map','.pathways','Why this framework exists','Nests, affected turtles, and conservation outcomes begin in different units. Follow these two pathways: monitoring estimates population status; demographic translation expresses impacts in a comparable currency.'),
    step('map','.currency','Meet in annual nesting females','Annual nesting females are the central population currency. ANE expresses expected contribution in that currency, allowing Take and Give schedules to be interpreted against the same population status.'),
    step('map','.gateway','Start where your data enter','Here we select incomplete monthly counts. Missing monitoring requires reconstruction before annual nests can be converted to annual nesters. An observed zero is different from no monitoring.',()=>{click('[data-source="nesting"]');click('[data-entry="0"]');}),
    step('map','.gateway','Different inputs join the same pathway','Now select annual nesting-female estimates. These enter after both monthly imputation and clutch-frequency conversion. The population model is still needed to estimate status, trend, and uncertainty.',()=>click('[data-entry="3"]')),
    step('monitoring','#stress-message','Missing months challenge reconstruction','The retained Bayesian Fourier nesting-season model reconstructs missing monthly monitoring. This graphic is schematic: it teaches the inference problem rather than showing a fitted population result.',()=>click('[data-stress="months"]')),
    step('monitoring','#stress-message','Skipped years weaken population-state recovery','Long gaps between monitoring years make the underlying population state harder to recover. This is a different problem from reconstructing months within a monitored nesting season.',()=>click('[data-stress="years"]')),
    step('monitoring','#stress-message','Short histories weaken trend precision','A short time series provides less information about trend. The final modeled annual-nester posterior starts projections; the separate RI/4 adult-female reporting quantity does not.',()=>click('[data-stress="short"]')),
    step('impact','#ane-chart','A small turtle has a delayed contribution','At 25% of asymptotic length (about 35.7 cm in this demonstration), survival and time to first nesting constrain the expected annual-nester contribution. Read the annual schedule and the cumulative total below it. These are fixed-parameter teaching calculations.',()=>set('size',25)),
    step('impact','#ane-chart','The same count can mean a different loss','Now the starting size is 99% of asymptotic length (about 141.3 cm), near first nesting in this demonstration. Expected ANE is larger and arrives earlier. First nesting is allocated once: the first-event probability includes the probability of not nesting previously.',()=>set('size',99)),
    step('give','#give-accounting','Conservation begins with a counterfactual','One protected nest × 77.9 eggs × (0.70 − 0.60 emergence) gives 7.79 additional emerged hatchlings. That incremental cohort must still pass through survival and first nesting before becoming Give ANE.',()=>click('[data-give="nest"]')),
    step('give','#give-accounting','Treatment alone is not benefit','Set counterfactual emergence to 0.70, equal to intervention emergence. Additional hatchlings and Give become zero. Give measures what would not otherwise have existed, rather than the total number treated.',()=>{click('[data-give="nest"]');set('counter',.7);}),
    step('give','#give-accounting','Biological entry point controls timing','Adult protection starts with deaths genuinely prevented. Its annual-nester contribution can be immediate or near-immediate; hatchling benefits arrive through a much longer demographic pathway.',()=>click('[data-give="adult"]')),
    step('projection','.steps','The order is scientific accounting','Subtract Take first, apply background growth and the retained abundance-scale √Q process error, then add Give and enforce non-negative abundance. Give has already been translated to the census year when it enters annual-nester currency.'),
    step('projection','#future-ledger','Give every scenario the same random future','Scenarios share the starting posterior draw, paired U/Q values, and the same future process realization z. The first transition uses the selected posterior row; later rows are sampled while preserving U/Q pairing. This ledger is a synthetic illustration.'),
    step('projection','#difference-chart','Compare paired differences','Read each scenario relative to status quo. Because the random future is shared, differences arise from Take and Give. Here the teaching schedules use Take = 3 and Give = 2 annual nesters per year.',()=>{set('take-amount',3);set('gain-amount',2);}),
    step('uncertainty','.uncertainty-list','Keep uncertainty and validation distinct','Population-state and trend uncertainty are retained. Baseline demographic parameters are fixed; changing demographic schedules are planned. Mathematical checks test implementation, not the biological appropriateness of your parameter choices.'),
    step('map','.currency','Use the map to explain; use Shiny to fit','Monitoring establishes population status. Demographic translation makes threats and incremental conservation benefits comparable to that status. Matched futures isolate their consequences. Finish to restore your original view and settings.')
  ];
  let index=0, saved=null, box=null, shade=null, ring=null;
  function capture() {
    return {hash:location.hash, panel:$('.panel:not([hidden])')?.id || 'map', x:window.scrollX,y:window.scrollY,
      source:$('[data-source][aria-pressed="true"]')?.dataset.source,
      entry:$('[data-entry][aria-pressed="true"]')?.dataset.entry,
      stress:$('[data-stress][aria-pressed="true"]')?.dataset.stress,
      give:$('[data-give][aria-pressed="true"]')?.dataset.give,
      inputs:all('.science input').map(e=>[e.id,e.value]),
      details:all('.science details').map(e=>[e,e.open]), focus:document.activeElement};
  }
  function close() {
    if(!saved)return;
    const s=saved; saved=null;
    box.remove();shade.remove();ring.remove();box=null;
    if(s.source)click(`[data-source="${s.source}"]`);
    if(s.entry!==undefined)click(`[data-entry="${s.entry}"]`);
    if(s.stress)click(`[data-stress="${s.stress}"]`);
    if(s.give)click(`[data-give="${s.give}"]`);
    s.inputs.forEach(([id,value])=>set(id,value));
    s.details.forEach(([e,open])=>e.open=open);
    click(`[data-panel="${s.panel}"]`);
    history.replaceState(null,'',location.pathname+location.search+s.hash);
    s.focus?.focus({preventScroll:true}); window.scrollTo(s.x,s.y);
  }
  function position() {
    if(!box)return;
    const target=$(steps[index].target); if(!target)return;
    const r=target.getBoundingClientRect(), vw=window.innerWidth,vh=window.innerHeight;
    Object.assign(ring.style,{left:Math.max(4,r.left-5)+'px',top:(r.top-5)+'px',width:Math.min(r.width+10,vw-8)+'px',height:(r.height+10)+'px'});
    const width=Math.min(380,vw-24);box.style.width=width+'px';
    const h=box.getBoundingClientRect().height;
    let left=Math.min(Math.max(12,r.left),vw-width-12),top=r.bottom+14;
    if(top+h>vh-12)top=r.top-h-14;
    if(top<12){left=vw-width-12;top=vh-h-12;}
    if(vw<640){left=12;top=vh-h-12;}
    Object.assign(box.style,{left:left+'px',top:Math.max(12,top)+'px'});
  }
  function show() {
    const s=steps[index]; click(`[data-panel="${s.panel}"]`);s.action?.();
    const target=$(s.target);if(!target){close();return;}
    box.innerHTML=`<div class="tour-count">Scientific walkthrough · ${index+1} / ${steps.length}</div><h2 id="tour-title">${s.title}</h2><p id="tour-description">${s.text}</p><div class="tour-controls"><button type="button" data-tour="exit">Exit tour</button><span></span><button type="button" data-tour="back" ${index===0?'disabled':''}>Back</button><button type="button" data-tour="next">${index===steps.length-1?'Finish':'Next →'}</button></div>`;
    target.scrollIntoView?.({behavior:'instant',block:'center'});
    position();box.querySelector('[data-tour="next"]').focus({preventScroll:true});
  }
  start.addEventListener('click',()=>{
    if(saved)return;
    saved=capture();index=0;
    shade=document.createElement('div');shade.className='tour-shade';shade.setAttribute('aria-hidden','true');
    ring=document.createElement('div');ring.className='tour-ring';ring.setAttribute('aria-hidden','true');
    box=document.createElement('section');box.className='tour-dialog';box.setAttribute('role','dialog');box.setAttribute('aria-modal','true');box.setAttribute('aria-labelledby','tour-title');box.setAttribute('aria-describedby','tour-description');
    box.addEventListener('click',e=>{const action=e.target.closest('[data-tour]')?.dataset.tour;if(action==='exit'||action==='next'&&index===steps.length-1)close();else if(action==='next'){index++;show();}else if(action==='back'&&index>0){index--;show();}});
    document.body.append(shade,ring,box);show();
  });
  document.addEventListener('keydown',e=>{
    if(!saved)return;
    if(e.key==='Escape'){e.preventDefault();close();}
    if(e.key==='Tab'){
      const buttons=[...box.querySelectorAll('button:not([disabled])')];const first=buttons[0],last=buttons.at(-1);
      if(e.shiftKey&&(document.activeElement===first||!box.contains(document.activeElement))){e.preventDefault();last.focus();}
      else if(!e.shiftKey&&(document.activeElement===last||!box.contains(document.activeElement))){e.preventDefault();first.focus();}
    }
  },true);
  window.addEventListener('resize',position);window.addEventListener('scroll',position,{passive:true});
}
