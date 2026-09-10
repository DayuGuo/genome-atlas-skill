const D = __DATA__; const BLOODV3 = __BLOOD__; const UI = __UI__; const FIND = __FIND__; const PGX = __PGX__; const PRS_EN = __PRS_EN__;
let LANG = 'zh';
try { LANG = localStorage.getItem('dayu-lang') || (navigator.language && !navigator.language.startsWith('zh') ? 'en' : 'zh'); } catch(e) {}
const F = {};
const t = k => ((UI[k] || [k,k])[LANG === 'zh' ? 0 : 1]).replace(/\{(\w+)(:\+)?\}/g, (m,v,plus)=>{const x=F[v]; if(x===undefined) return m; return (plus && x>0 ? '+' : '') + x;});
const zh = () => LANG === 'zh';
const NAME = () => zh() ? (D.name_zh || 'Sample') : (D.name_en || 'Sample');
const TITLE = () => zh() ? (D.title_zh || D.name_zh || '基因组图鉴') : (D.title_en || D.name_en || 'Genome Atlas');
const NS='http://www.w3.org/2000/svg';
const el=(p,tag,a)=>{const n=document.createElementNS(NS,tag);for(const k in a)n.setAttribute(k,a[k]);p.appendChild(n);return n};
const txt=(p,a,s)=>{const n=el(p,'text',a);n.textContent=s;return n};
const tip=(n,s)=>{const q=document.createElementNS(NS,'title');q.textContent=s;n.appendChild(q)};
const V=name=>getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const C=()=>({ink:V('--ink'),muted:V('--muted'),lab:V('--lab'),faint:V('--faint'),grid:V('--grid'),track:V('--track'),data:V('--data'),data2:V('--data2'),fd:V('--faintdata'),hero:V('--hero'),bg:V('--bg')});
const fmt=n=>n.toLocaleString('en-US');
const CH = {}; const NUM={'font-weight':800};
const foot=(s,c,w,y,txt_)=>txt(s,{x:w/2,y,'text-anchor':'middle','font-size':7,'font-weight':600,fill:c.faint,'letter-spacing':'.12em',class:'fade',style:'animation-delay:.9s'},txt_.toUpperCase());
const reveal=(id,fn)=>{CH[id]=fn;const n=document.getElementById(id);n.style.cursor='pointer';n.addEventListener('click',()=>{n.innerHTML='';fn(n,C())})};
const K=D.kpi;
F.depth=K.depth_auto; F.x=K.depth_x; F.y=K.depth_y; F.mt=fmt(K.depth_mt); F.callable=K.callable_gb; F.pass=fmt(K.pass_records); F.snv=fmt(K.snv); F.indel=fmt(K.indel); F.titv=K.titv;
F.cmp=fmt(D.chip.compared); F.nocall=D.chip.nocall; F.nonp=D.chip.nocall_nopass; F.indel_n=D.chip.in_deletion;
F.ng=fmt(D.anc.n_global); F.ne=fmt(D.anc.n_eas); F.knn=Object.entries(D.knn_eas).map(([k,v])=>`${k} ${v}`).join(' · '); F.knng=Object.entries(D.knn_global).map(([k,v])=>`${k} ${v}`).join(' · ');
F.near=Object.keys(D.near_eas)[0]; F.nearg=Object.entries(D.near_global).slice(0,2).map(([k,v])=>`${k} ${v}`).join(' ≈ ');
F.cvn=fmt(D.clinvar_total); F.cvd=D.clinvar_date; F.lof=D.lof.all; F.lofr=D.lof.rare; F.lofh=D.lof.rare_hom;
F.rho_i=D.prs_rho.imputed; F.rho_s=D.prs_rho.subset; F.sv=fmt(D.sv_total); F.rohn=D.roh_stats.n; F.rohmb=D.roh_stats.total_mb; F.rohmax=D.roh_stats.max_mb;
F.yterm=D.y_terminal; F.yformed=fmt(D.ypath[D.ypath.length-1].formed); F.mthg=D.mt.hg; F.mtq=D.mt.quality; F.mtn=D.mt.found.length; F.mtp=D.mt.private.length;

const renderAll=()=>{document.documentElement.setAttribute('lang',LANG==='zh'?'zh-Hans':'en');
  document.querySelectorAll('[data-i18n]').forEach(e=>{e.innerHTML=t(e.getAttribute('data-i18n'))});
  document.getElementById('langbtn').textContent=t('toggle');
  document.title = TITLE();
  for(const id in CH){const n=document.getElementById(id);n.innerHTML='';CH[id](n,C())}
  renderKpi(); renderFindings(); renderPgx(); renderHla(); renderMisc(); if(typeof deepRender==='function') deepRender();};
document.getElementById('langbtn').addEventListener('click',()=>{LANG=zh()?'en':'zh';try{localStorage.setItem('dayu-lang',LANG)}catch(e){};renderAll()});
if(window.matchMedia){window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change',renderAll)}
const prsName=n=>zh()?(PRS_EN[n]||n):n;

function renderKpi(){
  document.getElementById('kpi_depth').innerHTML=K.depth_auto+'<small>×</small>';
  document.getElementById('kpi_var').innerHTML=(K.pass_records/1e6).toFixed(2)+'<small>M</small>';
  document.getElementById('kpi_call').innerHTML=K.callable_gb+'<small>Gb</small>';
  document.getElementById('kpi_titv').textContent=K.titv.toFixed(2);
}

/* ── 00 per-chromosome depth + PASS variants ── */
reveal('chromdepth',(s,c)=>{
  const rows=D.chrom, x0=50, top=14, rh=12.2, W1=300, W2=260;
  txt(s,{x:x0,y:6,'font-size':7,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'平均深度 0–40×':'MEAN DEPTH 0–40×');
  txt(s,{x:x0+W1+120,y:6,'font-size':7,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'PASS 变异数（每 Mb）':'PASS VARIANTS PER MB');
  const maxd=Math.max(...rows.filter(r=>r.chrom!=='MT').map(r=>r.depth)), maxr=Math.max(...rows.filter(r=>r.len>0&&r.chrom!=='MT').map(r=>r.n_pass/r.len*1e6));
  rows.forEach((r,i)=>{const y=top+i*rh+6, sex=['X','Y','MT'].includes(r.chrom);
    txt(s,{x:x0-8,y:y+3,'text-anchor':'end','font-size':8,'font-weight':700,fill:c.ink,'letter-spacing':'.06em'},r.chrom);
    el(s,'line',{x1:x0,y1:y,x2:x0+W1,y2:y,stroke:c.grid,'stroke-width':.6});
    const d=r.chrom==='MT'?40:r.depth, bw=W1*Math.min(d,40)/40;
    const b=el(s,'rect',{x:x0,y:y-3.6,width:bw,height:7.2,rx:1.5,fill:sex?c.data2:c.data,class:'fade',style:`animation-delay:${i*.03}s`}); tip(b,`chr${r.chrom}: ${r.depth}×`);
    if(r.chrom==='MT') el(s,'rect',{x:x0+W1-4,y:y-3.6,width:2,height:7.2,fill:c.bg});
    txt(s,{x:x0+bw+5,y:y+3,'font-size':8,'font-weight':800,fill:c.ink},r.depth+'×');
    const px=x0+W1+120, ratio=r.len>0?r.n_pass/r.len*1e6:0, bw2=W2*Math.min(ratio,maxr)/maxr;
    el(s,'line',{x1:px,y1:y,x2:px+W2,y2:y,stroke:c.grid,'stroke-width':.6});
    const b2=el(s,'rect',{x:px,y:y-3.6,width:r.chrom==='MT'?0:bw2,height:7.2,rx:1.5,fill:sex?c.data2:c.data,class:'fade',style:`animation-delay:${.2+i*.03}s`}); tip(b2,`chr${r.chrom}: ${fmt(r.n_pass)} PASS variants`);
    txt(s,{x:px+(r.chrom==='MT'?0:bw2)+5,y:y+3,'font-size':7.5,'font-weight':800,fill:c.muted},r.chrom==='MT'?fmt(r.n_pass)+(zh()?' 个':''):fmt(r.n_pass));
  });
  [0,10,20,30,40].forEach(v=>{const x=x0+W1*v/40,yb=top+rows.length*rh+8;el(s,'line',{x1:x,y1:yb,x2:x,y2:yb+4,stroke:c.faint,'stroke-width':.8});txt(s,{x,y:yb+13,'text-anchor':'middle','font-size':7,fill:c.faint},v+'×')});
  foot(s,c,900,348,zh()?'X 约为常染色体一半、Y 更低：男性 · 线粒体 4,850× 超出刻度':'X about half of autosomes, Y lower: male · mitochondria 4,850× off scale');
});


/* ── 01 Y chain (YFull, from O-F438 down) ── */
reveal('ychain',(s,c)=>{
  const P=D.ypath, n=P.length, x0=40, x1=490, y=100;
  el(s,'line',{x1:x0,y1:y,x2:x1,y2:y,stroke:c.ink,'stroke-width':1,class:'draw'});
  P.forEach((p,i)=>{const x=x0+(x1-x0)*i/(n-1), r=4+Math.sqrt(p.der)*3, last=i===n-1;
    const cir=el(s,'circle',{cx:x,cy:y,r,fill:last?c.hero:c.data,stroke:c.bg,'stroke-width':1.5,class:'pop'}); cir.style.animationDelay=(i*110)+'ms'; tip(cir,`${p.snp} · ${p.der} derived / ${p.anc} ancestral · formed ${fmt(p.formed)} ybp`);
    txt(s,{x,y:y-r-10,'text-anchor':'middle','font-size':last?11:9.5,'font-weight':last?800:700,fill:c.ink,...NUM},p.snp.replace('O-',''));
    txt(s,{x,y:y+r+14,'text-anchor':'middle','font-size':7,'font-weight':600,fill:c.muted},p.der+(zh()?' 位点':' SNP'));
    txt(s,{x,y:y+r+25,'text-anchor':'middle','font-size':7,fill:c.faint},'~'+fmt(p.formed)+(zh()?' 年前':' ybp'));
  });
  txt(s,{x:x0,y:y+80,'font-size':8,'font-weight':600,fill:c.muted},zh()?'上游（ISOGG 树，第一轮已确认）：O-M122 → M324 → P201 → P164 → M134 → M117 → F8 → F438':'Upstream (ISOGG tree, round 1): O-M122 → M324 → P201 → P164 → M134 → M117 → F8 → F438');
  txt(s,{x:x0,y:y+96,'font-size':8,'font-weight':600,fill:c.muted},zh()?'每级兄弟支系全部祖先态（CTS12122 0/8、F316 0/16、Y277267 0/3、F15823 0/1）':'Sister branches all ancestral at every step (CTS12122 0/8, F316 0/16, Y277267 0/3, F15823 0/1)');
  txt(s,{x:x0,y:y+130,'font-size':24,'font-weight':700,fill:c.ink,...NUM},D.y_terminal);
  txt(s,{x:x0+150,y:y+130,'font-size':9,'font-weight':600,fill:c.muted},'YFull '+D.versions.yfull+(zh()?' · 叶节点':' · leaf branch'));
  txt(s,{x:x0,y:y+150,'font-size':8,'font-weight':600,fill:c.muted},(zh()?'23 个定义位点深度 9–26×，全部衍生态，零祖先态':'23 defining SNPs at 9–26× depth, all derived, none ancestral'));
  foot(s,c,520,296,zh()?'圆面积 = 该级衍生态位点数 · 最深色 = 终端支系 · 年代为 YFull 估计':'circle area = derived sites at that step · darkest = terminal · dates are YFull estimates');
});

/* ── 01 mt strip ── */
reveal('mtstrip',(s,c)=>{
  const L=16569, x0=24, x1=376, y=110;
  el(s,'line',{x1:x0,y1:y,x2:x1,y2:y,stroke:c.ink,'stroke-width':1});
  [0,4000,8000,12000,16569].forEach(p=>{const x=x0+(x1-x0)*p/L;el(s,'line',{x1:x,y1:y+4,x2:x,y2:y+8,stroke:c.faint,'stroke-width':.8});txt(s,{x,y:y+18,'text-anchor':'middle','font-size':7,fill:c.faint},p===16569?'16,569':fmt(p))});
  const hvs=x0+(x1-x0)*16024/L; el(s,'rect',{x:hvs,y:y-30,width:x1-hvs,height:26,fill:c.fd,opacity:.5}); txt(s,{x:x1,y:y-34,'text-anchor':'end','font-size':7,'font-weight':600,fill:c.muted},'HVS-I');
  el(s,'rect',{x:x0,y:y-30,width:(x1-x0)*576/L,height:26,fill:c.fd,opacity:.5}); txt(s,{x:x0,y:y-34,'font-size':7,'font-weight':600,fill:c.muted},'D-loop');
  const all=D.mt.found.map(v=>[v,false]).concat(D.mt.private.map(v=>[v,true])).sort((a,b)=>parseFloat(a[0])-parseFloat(b[0]));
  all.forEach(([v,isP],i)=>{const pos=parseFloat(v), x=x0+(x1-x0)*pos/L, ln=i%2===0?-22:22;
    el(s,'line',{x1:x,y1:y,x2:x,y2:y+ln*.55,stroke:c.faint,'stroke-width':.7});
    const d=el(s,'circle',{cx:x,cy:y+ln*.55,r:3.2,fill:isP?c.bg:c.data,stroke:isP?c.ink:c.data,'stroke-width':1,class:'pop'}); d.style.animationDelay=(i*20)+'ms'; tip(d,v+(isP?(zh()?' · 私有':' · private'):(zh()?' · 定义位点':' · defining')));
  });
  const hx=x0+(x1-x0)*16093/L; el(s,'path',{d:`M${hx-5} ${y+34} L${hx} ${y+27} L${hx+5} ${y+34} Z`,fill:c.hero}); txt(s,{x:hx,y:y+45,'text-anchor':'middle','font-size':7,'font-weight':700,fill:c.ink},'16093C 90%');
  txt(s,{x:x0,y:y+80,'font-size':22,'font-weight':700,fill:c.ink,...NUM},D.mt.hg);
  txt(s,{x:x0+80,y:y+80,'font-size':10,'font-weight':600,fill:c.muted},'PhyloTree '+D.versions.phylotree+' · '+D.mt.found.length+(zh()?' 个定义位点命中 · 质量 ':' defining sites · quality ')+D.mt.quality);
  const lines=zh()?['未命中的支系定义位点：'+D.mt.notfound.join('、'),'私有变异（'+D.mt.private.length+'）：'+D.mt.private.join(' '),'异质性：只有 m.16093T>C（90%，已知高变热点）；m.310/3107 为参考序列伪影']:
    ['Defining sites not found: '+D.mt.notfound.join(', '),'Private ('+D.mt.private.length+'): '+D.mt.private.join(' '),'Heteroplasmy: only m.16093T>C (90%, a known hypervariable site); m.310/3107 are reference artefacts'];
  lines.forEach((tx,i)=>txt(s,{x:x0,y:y+104+i*14,'font-size':7.8,'font-weight':500,fill:c.muted},tx));
  foot(s,c,400,296,zh()?'实心 = 支系定义位点 · 空心 = 私有变异 · 三角 = 异质性位点 · 4,850× 深度':'filled = branch-defining · hollow = private · triangle = heteroplasmic · 4,850× depth');
});

/* ── 02 scatter ── */
function scatter(s,c,pts,me,opts){
  const {x0=44,y0=14,W=360,H=290,xl='PC1',yl='PC2',color,labels}=opts;
  const xs=pts.map(p=>p.x).concat(me.x), ys=pts.map(p=>p.y).concat(me.y);
  let xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
  const dx=(xmax-xmin)*.06,dy=(ymax-ymin)*.06; xmin-=dx;xmax+=dx;ymin-=dy;ymax+=dy;
  const X=v=>x0+W*(v-xmin)/(xmax-xmin), Y=v=>y0+H*(1-(v-ymin)/(ymax-ymin));
  el(s,'rect',{x:x0,y:y0,width:W,height:H,fill:'none',stroke:c.grid,'stroke-width':.6});
  txt(s,{x:x0+W/2,y:y0+H+22,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted,'letter-spacing':'.1em'},xl);
  txt(s,{x:x0-30,y:y0+H/2,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted,'letter-spacing':'.1em',transform:`rotate(-90 ${x0-30} ${y0+H/2})`},yl);
  pts.forEach(p=>{const d=el(s,'circle',{cx:X(p.x).toFixed(1),cy:Y(p.y).toFixed(1),r:2.1,fill:color(p),opacity:.85});tip(d,p.t)});
  if(labels){const groups={};pts.forEach(p=>{(groups[p.g]=groups[p.g]||[]).push(p)});
    Object.entries(groups).forEach(([g,arr])=>{const mx=arr.reduce((a,p)=>a+p.x,0)/arr.length,my=arr.reduce((a,p)=>a+p.y,0)/arr.length;
      txt(s,{x:X(mx),y:Y(my)-6,'text-anchor':'middle','font-size':7.5,'font-weight':700,fill:c.ink,stroke:c.bg,'stroke-width':2.5,'paint-order':'stroke','letter-spacing':'.06em'},g.toUpperCase())})}
  const mx=X(me.x),my=Y(me.y);
  el(s,'rect',{x:mx-5.5,y:my-5.5,width:11,height:11,fill:c.hero,stroke:c.bg,'stroke-width':1.5,transform:`rotate(45 ${mx} ${my})`,class:'pop'});
  txt(s,{x:mx+9,y:my+3,'font-size':9,'font-weight':800,fill:c.ink,stroke:c.bg,'stroke-width':3,'paint-order':'stroke'},NAME());
  if(opts.foot) foot(s,c,420,336,opts.foot);
}
reveal('pcaglobal',(s,c)=>{const pts=D.pca_global.pts.map(p=>({g:p[0],x:p[2],y:p[3],t:`${p[1]} (${p[0]})`}));
  scatter(s,c,pts,{x:D.pca_global.me[0],y:D.pca_global.me[1]},{color:p=>p.g==='EAS'?c.data:c.fd,labels:true,foot:zh()?'一点 = 一个人 · 深蓝 = 东亚 · 菱形 = '+NAME():'one dot = one person · dark blue = East Asian · diamond = '+NAME()})});
reveal('pcaeas',(s,c)=>{const shade={CHB:c.hero,CHS:c.data,JPT:c.data2,CDX:c.fd,KHV:c.fd};const pts=D.pca_eas.pts.map(p=>({g:p[0],x:p[1],y:p[2],t:p[0]}));
  scatter(s,c,pts,{x:D.pca_eas.me[0],y:D.pca_eas.me[1]},{color:p=>shade[p.g]||c.fd,labels:true,foot:zh()?'一点 = 一个人 · 明度 = 人群 · 菱形 = '+NAME():'one dot = one person · shade = population · diamond = '+NAME()})});

/* ── 03 findings ── */
function renderFindings(){
  document.getElementById('findings').innerHTML=FIND.map(f=>`<div class="finding"><div><div class="g">${f.gene}<small>${zh()?f.gt_zh:f.gt_en}</small></div><span class="stamp${f.hot?' hot':''}">${zh()?f.tag_zh:f.tag_en}</span></div><div><div class="b">${zh()?f.zh:f.en}</div>${f.act_zh?`<div class="act">${zh()?f.act_zh:f.act_en}</div>`:''}</div></div>`).join('');
}
/* ── 03 ClinVar bars ── */
reveal('cvbars',(s,c)=>{
  const order=['Benign','Likely_benign','Uncertain_significance','Conflicting_classifications_of_pathogenicity','drug_response','association','risk_factor','other','not_provided','protective','Pathogenic/Likely_pathogenic','Pathogenic','Likely_pathogenic'];
  const rows=order.filter(k=>D.clinvar[k]).map(k=>[k,D.clinvar[k]]); const x0=250,x1=880,top=10,rh=13, max=rows[0][1];
  const nm={Benign:['良性','Benign'],Likely_benign:['可能良性','Likely benign'],Uncertain_significance:['意义未明','Uncertain'],Conflicting_classifications_of_pathogenicity:['分类冲突','Conflicting'],drug_response:['药物反应','Drug response'],association:['关联','Association'],risk_factor:['风险因素','Risk factor'],other:['其它','Other'],not_provided:['未提供','Not provided'],protective:['保护性','Protective'],'Pathogenic/Likely_pathogenic':['致病/可能致病','P/LP'],Pathogenic:['致病','Pathogenic'],Likely_pathogenic:['可能致病','Likely pathogenic']};
  rows.forEach(([k,v],i)=>{const y=top+i*rh+5, w=Math.max(1.5,(x1-x0)*Math.sqrt(v)/Math.sqrt(max)), pat=/athogenic/.test(k)&&!/Conflicting/.test(k);
    txt(s,{x:x0-8,y:y+3,'text-anchor':'end','font-size':8,'font-weight':pat?800:600,fill:pat?c.ink:c.muted},(nm[k]||[k,k])[zh()?0:1]);
    const b=el(s,'rect',{x:x0,y:y-3.5,width:w,height:7,rx:1.5,fill:pat?c.hero:(/Conflicting|Uncertain/.test(k)?c.data:c.data2),class:'fade',style:`animation-delay:${i*.04}s`}); tip(b,`${k}: ${fmt(v)}`);
    txt(s,{x:x0+w+5,y:y+3,'font-size':8,'font-weight':800,fill:c.ink},fmt(v));
  });
  foot(s,c,900,186,zh()?'条长按平方根压缩 · 最深色 = 致病/可能致病 · 只统计 PASS 变异':'bar length on a square-root scale · darkest = pathogenic / likely pathogenic · PASS variants only');
});

/* ── 04 PGx ── */
function renderPgx(){
  document.getElementById('pgx').innerHTML=PGX.map(r=>`<tr><td class="gene">${r.gene}</td><td class="g">${r.dip}</td><td>${zh()?r.ph_zh:r.ph_en}</td><td>${zh()?r.zh:r.en}</td><td class="rs">${r.src}</td></tr>`).join('');
}
function renderHla(){
  const order=['HLA-A','HLA-B','HLA-C','HLA-DRB1','HLA-DQA1','HLA-DQB1','HLA-DPA1','HLA-DPB1'];
  const hot=new Set(['B*13:01:01']);
  document.getElementById('hla').innerHTML=order.map(g=>{const v=D.hla[g]; if(!v) return ''; const f=a=>hot.has(a)?`<em>${a.replace(/:\d+$/,'')}</em>`:a.replace(/(\*\d+:\d+).*/,'$1'); return `<div><div class="k">${g}</div><div class="v">${f(v[0])} / ${f(v[1])}</div></div>`}).join('');
}

/* ── 05 PRS three versions ── */
reveal('prs',(s,c)=>{
  const rows=[...D.prs].sort((a,b)=>b.pct_EAS-a.pct_EAS), x0=250, x1=850, top=16, rh=21, X=v=>x0+(x1-x0)*v/100;
  el(s,'rect',{x:X(25),y:top-8,width:X(75)-X(25),height:rows.length*rh+4,fill:c.track,opacity:.6});
  [0,25,50,75,100].forEach(v=>{el(s,'line',{x1:X(v),y1:top-8,x2:X(v),y2:top+rows.length*rh-2,stroke:v===50?c.faint:c.grid,'stroke-width':v===50?.9:.5});txt(s,{x:X(v),y:top+rows.length*rh+12,'text-anchor':'middle','font-size':7.5,fill:c.faint},v)});
  rows.forEach((r,i)=>{const y=top+i*rh+6, hi=r.pct_EAS>=90||r.pct_EAS<=10;
    txt(s,{x:x0-12,y:y+3,'text-anchor':'end','font-size':9,'font-weight':hi?800:600,fill:hi?c.ink:c.lab},prsName(r.trait));
    el(s,'line',{x1:X(r.pct_Han),y1:y,x2:X(r.pct_EAS),y2:y,stroke:c.data,'stroke-width':1,opacity:.5}); const q=el(s,'circle',{cx:X(r.pct_Han),cy:y,r:3.6,fill:c.bg,stroke:c.data,'stroke-width':1.3});tip(q,(zh()?'汉族 208 人中 ':'among 208 Han ')+r.pct_Han)
    const mx=X(r.pct_EAS); const d=el(s,'rect',{x:mx-4.5,y:y-4.5,width:9,height:9,fill:c.hero,stroke:c.bg,'stroke-width':1.2,transform:`rotate(45 ${mx} ${y})`,class:'pop'}); d.style.animationDelay=(i*25)+'ms'; tip(d,`${r.trait}: WGS ${r.pct_EAS} (Han ${r.pct_Han}) · coverage ${r.coverage_pct}% · z ${r.z_vs_EAS}`);
    txt(s,{x:x1+10,y:y+3,'font-size':8.5,'font-weight':800,fill:c.ink},r.pct_EAS.toFixed(0));
    txt(s,{x:x1+34,y:y+3,'font-size':7.5,fill:c.faint},r.coverage_pct+'%');
  });
  txt(s,{x:x1+10,y:top-10,'font-size':7,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'百分位 · 覆盖':'PCT · COVERAGE');
  foot(s,c,900,696,zh()?'菱形 = 在 504 名东亚人中的百分位 · 空心圆 = 在 208 名汉族中 · 灰带 = 中间一半的人':'diamond = percentile among 504 East Asians · hollow = among 208 Han · band = middle half');
});

/* ── 06 CNV depth panel ── */
reveal('cnv',(s,c)=>{
  const rows=D.cnv, x0=120, x1=470, top=14, rh=24, X=v=>x0+(x1-x0)*Math.min(v,2.2)/2.2;
  [0,0.5,1,1.5,2].forEach(v=>{el(s,'line',{x1:X(v),y1:top-6,x2:X(v),y2:top+rows.length*rh,stroke:v===1?c.faint:c.grid,'stroke-width':v===1?1:.5,'stroke-dasharray':v===1?'':'2 3'});txt(s,{x:X(v),y:top+rows.length*rh+12,'text-anchor':'middle','font-size':7.5,fill:c.faint},v===1?(zh()?'1.0 = 两份':'1.0 = two copies'):v.toFixed(1))});
  rows.forEach((r,i)=>{const y=top+i*rh+8, w=X(r.ratio)-x0;
    txt(s,{x:x0-8,y:y+3,'text-anchor':'end','font-size':9,'font-weight':r.copies<2?800:600,fill:c.ink},r.locus);
    const b=el(s,'rect',{x:x0,y:y-5,width:Math.max(w,1.5),height:10,rx:2,fill:r.copies===0?c.hero:(r.copies===1?c.data:c.data2),class:'fade',style:`animation-delay:${i*.05}s`}); tip(b,`${r.locus}: depth ratio ${r.ratio} · ${r.copies} copies · ${r.evidence}`);
    txt(s,{x:X(r.ratio)+6,y:y+3,'font-size':8,'font-weight':800,fill:c.ink},r.copies+(zh()?' 份':(r.copies===1?' copy':' copies')));
  });
  foot(s,c,520,326,zh()?'条长 = 该区读段深度 / 全基因组平均 · 最深色 = 零份（纯合缺失）· SMN 来自专用工具':'bar = read depth in the locus / genome average · darkest = zero copies (homozygous deletion) · SMN from a dedicated caller');
});
/* ── 06 STR panel ── */
reveal('str',(s,c)=>{
  const rows=D.str.filter(r=>r.thr).sort((a,b)=>b.max/b.thr-a.max/a.thr).slice(0,22), x0=78, x1=372, top=14, rh=13.6, X=v=>x0+(x1-x0)*Math.min(v,1.1)/1.1;
  el(s,'line',{x1:X(1),y1:top-6,x2:X(1),y2:top+rows.length*rh,stroke:c.ink,'stroke-width':1}); txt(s,{x:X(1),y:top-9,'text-anchor':'middle','font-size':7,'font-weight':700,fill:c.ink},zh()?'致病阈值':'PATHOGENIC');
  rows.forEach((r,i)=>{const y=top+i*rh+6, w=X(r.max/r.thr)-x0;
    txt(s,{x:x0-6,y:y+3,'text-anchor':'end','font-size':8,'font-weight':700,fill:c.ink},r.locus);
    el(s,'line',{x1:x0,y1:y,x2:X(1),y2:y,stroke:c.grid,'stroke-width':.5});
    const b=el(s,'rect',{x:x0,y:y-3.5,width:Math.max(w,1),height:7,rx:1.5,fill:c.data,class:'fade',style:`animation-delay:${i*.03}s`}); tip(b,`${r.locus} (${r.unit}): ${r.gt} repeats · threshold ≥${r.thr}`);
    txt(s,{x:x0+w+4,y:y+3,'font-size':7.5,'font-weight':800,fill:c.ink},r.gt);
  });
  foot(s,c,400,326,zh()?'条长 = 最长等位基因 / 致病阈值 · 全部在阈值左侧 · ExpansionHunter 5':'bar = longest allele / pathogenic threshold · all left of the line · ExpansionHunter 5');
});

/* ── 07 ROH ideogram ── */
reveal('roh',(s,c)=>{
  const chr=[...Array(22)].map((_,i)=>String(i+1)), x0=40, top=12, rh=12.4, maxL=D.chrlen['1'], W=830;
  chr.forEach((ch,i)=>{const y=top+i*rh+5, L=D.chrlen[ch], w=W*L/maxL;
    txt(s,{x:x0-8,y:y+3,'text-anchor':'end','font-size':8,'font-weight':700,fill:c.ink},ch);
    el(s,'rect',{x:x0,y:y-2.5,width:w,height:5,rx:2.5,fill:c.track});
    const cen=D.cen[ch]; if(cen) el(s,'line',{x1:x0+W*cen*1e6/maxL,y1:y-5,x2:x0+W*cen*1e6/maxL,y2:y+5,stroke:c.faint,'stroke-width':.9});
    D.roh.filter(r=>r.chrom===ch).forEach(r=>{const rx=x0+W*r.start/maxL, rw=Math.max(2,W*(r.end-r.start)/maxL); const b=el(s,'rect',{x:rx,y:y-3.2,width:rw,height:6.4,rx:1,fill:r.mb>=1.5?c.hero:c.data,class:'pop'}); tip(b,`chr${ch}:${fmt(r.start)}-${fmt(r.end)} · ${r.mb} Mb · q ${r.q}`)});
  });
  foot(s,c,900,296,zh()?'深色段 = ≥1.5 Mb 的纯合区 · 浅色 = 1–1.5 Mb · 小刻度 = 着丝粒 · 跨着丝粒/端粒的伪影段已剔除':'dark = homozygous run ≥1.5 Mb · light = 1–1.5 Mb · tick = centromere · segments spanning centromeres/telomeres removed');
});
function renderMisc(){
  const items=zh()?[['SMN1 / SMN2',`${D.smn.SMN1} / ${D.smn.SMN2} 份 · 非 SMA 携带者`],['线粒体异质性','m.16093T>C 90%'],['结构变异总数',`${fmt(D.sv_total)}（DEL ${fmt(D.sv_counts.DEL)} · DUP ${D.sv_counts.DUP} · INV ${D.sv_counts.INV}）`],['罕见功能丧失变异',`${D.lof.rare} 个 · 纯合 ${D.lof.rare_hom} 个均在重复区`],['X 染色体','PAR 重调后 6,919 个杂合位点'],['近亲检查',`ROH ≥1 Mb ${D.roh_stats.n} 段 · 合计 ${D.roh_stats.total_mb} Mb · 无 >5 Mb`]]:
    [['SMN1 / SMN2',`${D.smn.SMN1} / ${D.smn.SMN2} copies · not an SMA carrier`],['mt heteroplasmy','m.16093T>C 90%'],['Structural variants',`${fmt(D.sv_total)} (DEL ${fmt(D.sv_counts.DEL)} · DUP ${D.sv_counts.DUP} · INV ${D.sv_counts.INV})`],['Rare loss-of-function',`${D.lof.rare} · ${D.lof.rare_hom} homozygous, all in repeats`],['X chromosome','6,919 heterozygous PAR sites after re-calling'],['Relatedness',`${D.roh_stats.n} ROH ≥1 Mb · ${D.roh_stats.total_mb} Mb total · none >5 Mb`]];
  document.getElementById('misc').innerHTML=items.map(([k,v])=>`<div><div class="k">${k}</div><div class="v" style="font-size:12px">${v}</div></div>`).join('');
}

/* ══════════ deep-dive additions (phase 6) ══════════ */
const AC=()=>({n:V('--ancn'),s:V('--ancs'),a:V('--arch')});
const CHR22=[...Array(22)].map((_,i)=>String(i+1));
const calCHB=D.la_calib.find(x=>x.pop==='CHB'), calCHS=D.la_calib.find(x=>x.pop==='CHS');
const CALCH=['1','2','6','22'];
const calLen=CALCH.reduce((a,c)=>a+D.chrlen[c],0);
const dayuCal=CALCH.reduce((a,c)=>a+D.la_per_chrom.find(x=>String(x.chrom)===c).NorthEA*D.chrlen[c],0)/calLen;
F.north=(D.la_global.NorthEA*100).toFixed(1);
F.south=(D.la_global.SouthEA*100).toFixed(1);
F.noise=((D.la_global.European+D.la_global.SouthAsian)*100).toFixed(1);
F.chb=(calCHB.north_mean*100).toFixed(0); F.chs=(calCHS.north_mean*100).toFixed(0);
F.dayucal=(dayuCal*100).toFixed(1);
F.sdchb=((dayuCal-calCHB.north_mean)/calCHB.north_sd).toFixed(1);
F.archmb=D.archaic_summary.span_mb; F.neamb=D.archaic_summary.neanderthal_mb; F.denmb=D.archaic_summary.denisovan_mb;
F.archn=D.archaic_summary.merged; F.archpct=(D.archaic_summary.span_mb/2875*100).toFixed(1);
F.phhet=fmt(D.phase.het); F.phased=fmt(D.phase.phased); F.phpct=(D.phase.phased/D.phase.het*100).toFixed(1);
F.n50=(D.phase.n50_kb/1000).toFixed(1); F.blocks=fmt(D.phase.blocks);
F.mtcn=D.somatic.mtDNA_copies_per_cell; F.ydr=D.somatic.Y_depth_ratio; F.telk7=fmt(D.telomere.k7);
F.snvn=fmt(D.spectrum.reduce((a,b)=>a+b.n,0));
F.cpg=(D.spectrum.filter(x=>/\[C>T\]G/.test(x.ctx)).reduce((a,b)=>a+b.frac,0)*100).toFixed(1);

const ANC_ZH={"China_Henan_Jiaozuoniecunsite_LBA_IA":"河南 焦作聂村 晚商–铁器时代","China_Shandong_Chengziya_Yueshi":"山东 城子崖 岳石文化","China_Henan_Pingliangtaisite_LN":"河南 平粮台 龙山文化","China_Baligang_BA_EasternZhou":"河南 八里岗 东周","China_Jinan_LiuJiaZhuang_Shang":"山东 济南刘家庄 商","China_Jining_YinJiaCheng_Longshan":"山东 济宁尹家城 龙山","China_Shandong_Chengziya_Longshan":"山东 城子崖 龙山","China_InnerMongolia_Erdaojingzi_LN":"内蒙古 二道井子 新石器晚期","China_Henan_Haojiatai_LN":"河南 郝家台 龙山","China_Baligang_LN_Longshan":"河南 八里岗 龙山","China_Shandong_Dinggong_LN":"山东 丁公 龙山","China_Henan_Wadiansite_LN":"河南 瓦店 龙山","China_Baligang_LN_Shijiahe":"河南 八里岗 石家河","China_Qingdao_BeiQian_Dawenkou":"山东 青岛北阡 大汶口","China_Baligang_LN_Yangshao":"河南 八里岗 仰韶","China_Shanxi_Shengedaliang_LN":"陕西 神圪垯梁 龙山","China_LBA_EIA":"中国 晚青铜–早铁器","China_Qinghai_Dacaozisite_IA":"青海 大草子 铁器时代","China_Qinghai_Lajiasite_LN":"青海 喇家 龙山","Japan_KofunPeriod":"日本 古坟时代"};

/* ── chromosome painting ── */
reveal('painting',(s,c)=>{
  const a=AC(), x0=34, rows=11, rh=44, maxL=D.chrlen['1'], COLW=390;
  const seg={}; D.la_segments.forEach(g=>{(seg[g.chrom]=seg[g.chrom]||[]).push(g)});
  CHR22.forEach((ch,i)=>{
    const col=i<rows?0:1, row=i%rows, bx=col===0?x0:x0+470, y=18+row*rh, w=COLW*D.chrlen[ch]/maxL;
    txt(s,{x:bx-6,y:y+11,'text-anchor':'end','font-size':10,'font-weight':700,fill:c.ink},ch);
    [0,1].forEach(h=>{const yy=y+h*11;
      el(s,'rect',{x:bx,y:yy,width:w,height:9,rx:1.5,fill:a.n,opacity:.5});
      (seg[ch]||[]).filter(g=>g.hap===h+1).forEach(g=>{
        const gx=bx+w*g.start/D.chrlen[ch], gw=Math.max(1.2,w*(g.end-g.start)/D.chrlen[ch]);
        const r=el(s,'rect',{x:gx,y:yy,width:gw,height:9,rx:1,fill:g.anc==='SouthEA'?a.s:c.fd});
        tip(r,`chr${ch}:${fmt(g.start)}-${fmt(g.end)} · ${g.anc}`)})});
    const cen=D.cen[ch]; if(cen) el(s,'line',{x1:bx+w*cen*1e6/D.chrlen[ch],y1:y-2,x2:bx+w*cen*1e6/D.chrlen[ch],y2:y+22,stroke:c.bg,'stroke-width':1.6});
    const so=D.la_per_chrom.find(r=>String(r.chrom)===ch);
    if(so) txt(s,{x:bx+w+7,y:y+13,'font-size':9,'font-weight':700,fill:c.muted},(so.SouthEA*100).toFixed(0)+'%');
  });
  foot(s,c,900,496,zh()?'每条染色体两行 = 父母各给的一条单倍型 · 只画 ≥0.5 Mb 的非北方片段 · 右侧 = 该染色体的南方比例':'two rows per chromosome = the two parental haplotypes · only non-northern segments ≥0.5 Mb · right = southern share');
});

/* ── calibration ── */
reveal('calib',(s,c)=>{
  const a=AC(), x0=170,x1=770,y0=34,rh=32;
  const rows=[{lab:NAME(),n:dayuCal,sd:0,me:true},
              {lab:zh()?'北京汉 CHB · 20 人':'Beijing Han (CHB), n=20',n:calCHB.north_mean,sd:calCHB.north_sd},
              {lab:zh()?'南方汉 CHS · 20 人':'southern Han (CHS), n=20',n:calCHS.north_mean,sd:calCHS.north_sd}];
  const X=v=>x0+(x1-x0)*(v-0.5)/0.5;
  [0.5,0.6,0.7,0.8,0.9,1].forEach(v=>{el(s,'line',{x1:X(v),y1:y0-12,x2:X(v),y2:y0+rows.length*rh-10,stroke:c.grid,'stroke-width':.7});
    txt(s,{x:X(v),y:y0+rows.length*rh+6,'text-anchor':'middle','font-size':8,fill:c.faint},(v*100)+'%')});
  txt(s,{x:x0,y:14,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'1、2、6、22 号染色体的北方东亚成分':'NORTHERN EAST ASIAN SHARE, CHROMOSOMES 1, 2, 6, 22');
  rows.forEach((r,i)=>{const y=y0+i*rh;
    txt(s,{x:x0-12,y:y+4,'text-anchor':'end','font-size':10.5,'font-weight':r.me?800:600,fill:r.me?c.ink:c.lab},r.lab);
    if(r.sd>0){el(s,'line',{x1:X(r.n-r.sd),y1:y,x2:X(r.n+r.sd),y2:y,stroke:c.faint,'stroke-width':1.4});
      [-1,1].forEach(k=>el(s,'line',{x1:X(r.n+k*r.sd),y1:y-4,x2:X(r.n+k*r.sd),y2:y+4,stroke:c.faint,'stroke-width':1.4}))}
    const d=el(s,'circle',{cx:X(r.n),cy:y,r:r.me?6:4.5,fill:r.me?a.n:c.bg,stroke:r.me?c.bg:c.lab,'stroke-width':r.me?2:1.5,class:'pop'});
    tip(d,`${(r.n*100).toFixed(1)}%${r.sd?` ± ${(r.sd*100).toFixed(1)}`:''}`);
    txt(s,{x:X(r.n+(r.sd||0))+15,y:y+4,'font-size':10.5,'font-weight':800,fill:c.ink},(r.n*100).toFixed(1)+'%')});
  foot(s,c,900,166,zh()?'横线 = 组内 ±1 标准差 · 参考个体打分时已从参考面板移出':'whisker = ±1 sd within the group · reference individuals held out of the panel when scored');
});

/* ── Human Origins PCA with ancient projection ── */
reveal('hopca',(s,c)=>{
  const a=AC(), x0=44,y0=16,W=440,H=390;
  const pts=D.ho_modern, anc=D.ho_ancient_pts, me=D.ho_me;
  const xs=pts.map(p=>p[1]).concat(anc.map(p=>p.pc1),[me[0]]), ys=pts.map(p=>p[2]).concat(anc.map(p=>p.pc2),[me[1]]);
  let xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
  const dx=(xmax-xmin)*.05,dy=(ymax-ymin)*.05; xmin-=dx;xmax+=dx;ymin-=dy;ymax+=dy;
  const X=v=>x0+W*(v-xmin)/(xmax-xmin), Y=v=>y0+H*(1-(v-ymin)/(ymax-ymin));
  el(s,'rect',{x:x0,y:y0,width:W,height:H,fill:'none',stroke:c.grid,'stroke-width':.6});
  const COL={Han:a.n,SouthEA:a.s,NorthAsia:c.data2};
  pts.forEach(p=>{el(s,'circle',{cx:X(p[1]).toFixed(1),cy:Y(p[2]).toFixed(1),r:2,fill:COL[p[0]]||c.fd,opacity:COL[p[0]]?.6:.34})});
  anc.forEach(p=>{const d=el(s,'rect',{x:X(p.pc1)-2.4,y:Y(p.pc2)-2.4,width:4.8,height:4.8,fill:'none',stroke:c.ink,'stroke-width':1,opacity:.75});
    tip(d,`${p.label} · ${fmt(p.date)} BP`)});
  const SHOW={Henan:[0,-26],Shandong:[46,-12],Guangdong:[-44,16],Fujian:[40,10],Sichuan:[-52,-4]};
  D.ho_prov.forEach(p=>{const off=SHOW[p.label]; if(!off) return;
    const lab=zh()?({Shandong:'山东',Henan:'河南',Fujian:'福建',Guangdong:'广东',Sichuan:'四川'}[p.label]):p.label;
    const px=X(p.pc1),py=Y(p.pc2);
    el(s,'circle',{cx:px,cy:py,r:2.6,fill:'none',stroke:c.ink,'stroke-width':1.2});
    el(s,'line',{x1:px,y1:py,x2:px+off[0],y2:py+off[1],stroke:c.faint,'stroke-width':.7});
    txt(s,{x:px+off[0]+(off[0]>=0?3:-3),y:py+off[1]+3,'text-anchor':off[0]>=0?'start':'end','font-size':9,'font-weight':700,fill:c.ink,stroke:c.bg,'stroke-width':2.6,'paint-order':'stroke'},lab)});
  const mx=X(me[0]),my=Y(me[1]);
  el(s,'rect',{x:mx-5.5,y:my-5.5,width:11,height:11,fill:c.hero,stroke:c.bg,'stroke-width':1.5,transform:`rotate(45 ${mx} ${my})`,class:'pop'});
  txt(s,{x:mx+11,y:my+4,'font-size':10,'font-weight':800,fill:c.ink,stroke:c.bg,'stroke-width':3,'paint-order':'stroke'},NAME());
  txt(s,{x:x0+W/2,y:y0+H+22,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted,'letter-spacing':'.1em'},'PC1');
  txt(s,{x:x0-30,y:y0+H/2,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted,'letter-spacing':'.1em',transform:`rotate(-90 ${x0-30} ${y0+H/2})`},'PC2');
  foot(s,c,520,462,zh()?'空心方块 = 古代个体（投影）· 圆点 = 现代个体 · 菱形 = '+NAME():'open squares = ancient (projected) · dots = present-day · diamond = '+NAME());
});

/* ── distance vs time ── */
reveal('timedist',(s,c)=>{
  const a=AC(), x0=54,y0=16,W=330,H=390;
  const rows=D.ho_near_ancient.filter(r=>r.date>0).slice(0,14);
  const dmax=Math.max(...rows.map(r=>r.d))*1.1, tmax=7600;
  const X=v=>x0+W*(1-v/tmax), Y=v=>y0+H*(v/dmax);
  el(s,'rect',{x:x0,y:y0,width:W,height:H,fill:'none',stroke:c.grid,'stroke-width':.6});
  [7000,5000,3000,1000].forEach(v=>{el(s,'line',{x1:X(v),y1:y0,x2:X(v),y2:y0+H,stroke:c.grid,'stroke-width':.5});
    txt(s,{x:X(v),y:y0+H+16,'text-anchor':'middle','font-size':8,fill:c.faint},fmt(v))});
  [0,0.01,0.02].forEach(v=>{if(v>dmax)return; el(s,'line',{x1:x0,y1:Y(v),x2:x0+W,y2:Y(v),stroke:c.grid,'stroke-width':.5});
    txt(s,{x:x0-8,y:Y(v)+3,'text-anchor':'end','font-size':8,fill:c.faint},v.toFixed(2))});
  const xs=rows.map(r=>r.date), ys=rows.map(r=>r.d);
  const mx=xs.reduce((p,q)=>p+q,0)/xs.length, my=ys.reduce((p,q)=>p+q,0)/ys.length;
  const b1=xs.reduce((p,x,i)=>p+(x-mx)*(ys[i]-my),0)/xs.reduce((p,x)=>p+(x-mx)**2,0), b0=my-b1*mx;
  el(s,'line',{x1:X(500),y1:Y(b0+b1*500),x2:X(7200),y2:Y(b0+b1*7200),stroke:c.faint,'stroke-width':1.2,'stroke-dasharray':'4 4'});
  rows.forEach(r=>{const d=el(s,'circle',{cx:X(r.date),cy:Y(r.d),r:3.6+Math.sqrt(r.n),fill:a.n,opacity:.78,stroke:c.bg,'stroke-width':1.2,class:'pop'});
    tip(d,`${r.label} · n=${r.n} · ${fmt(Math.round(r.date))} BP · d=${r.d.toFixed(4)}`)});
  const near=rows[0];
  txt(s,{x:X(near.date)+10,y:Y(near.d)+3,'font-size':9,'font-weight':700,fill:c.ink,stroke:c.bg,'stroke-width':2.6,'paint-order':'stroke'},zh()?'河南 焦作':'Henan Jiaozuo');
  txt(s,{x:x0+W/2,y:y0+H+34,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted,'letter-spacing':'.08em'},zh()?'年代（距今年数）':'YEARS BEFORE PRESENT');
  txt(s,{x:x0-38,y:y0+H/2,'text-anchor':'middle','font-size':8,'font-weight':700,fill:c.muted,'letter-spacing':'.08em',transform:`rotate(-90 ${x0-38} ${y0+H/2})`},zh()?'到样本的距离':'DISTANCE TO SAMPLE');
  foot(s,c,420,462,zh()?'圆越大样本越多 · 虚线为线性趋势':'bigger circle = more individuals · dashed line = linear trend');
});

/* ── archaic ideogram ── */
reveal('archaic',(s,c)=>{
  const x0=30,x1=880,rows=11,rh=28,maxL=D.chrlen['1'];
  const by={}; D.archaic.forEach(g=>{(by[g.chrom]=by[g.chrom]||[]).push(g)});
  CHR22.forEach((ch,i)=>{
    const col=i<rows?0:1, row=i%rows, bx=col===0?x0:x0+440, bw=(x1-x0-60)/2;
    const y=20+row*rh, w=bw*D.chrlen[ch]/maxL;
    txt(s,{x:bx-6,y:y+7,'text-anchor':'end','font-size':10,'font-weight':700,fill:c.ink},ch);
    el(s,'rect',{x:bx,y:y,width:w,height:10,rx:2,fill:c.track});
    const cen=D.cen[ch]; if(cen) el(s,'line',{x1:bx+w*cen*1e6/D.chrlen[ch],y1:y-2,x2:bx+w*cen*1e6/D.chrlen[ch],y2:y+12,stroke:c.faint,'stroke-width':1});
    (by[ch]||[]).forEach(g=>{
      const gx=bx+w*g.start/D.chrlen[ch], gw=Math.max(1.5,w*(g.end-g.start)/D.chrlen[ch]);
      const col2=g.src==='Denisovan'?c.data2:AC().a;
      const r=el(s,'rect',{x:gx,y:g.hom?y-1.5:y,width:gw,height:g.hom?13:10,rx:1.5,fill:col2,opacity:g.hom?1:.72});
      tip(r,`chr${ch}:${fmt(g.start)}-${fmt(g.end)} · ${g.src} · ${g.mb} Mb${g.hom?' · homozygous':''}`)});
    txt(s,{x:bx+w+6,y:y+8,'font-size':9,fill:c.muted},(by[ch]||[]).length);
  });
  foot(s,c,900,336,zh()?'绿色 = 尼安德特 · 浅蓝 = 丹尼索瓦 · 加高 = 两条染色体都带 · 右侧 = 片段数':'green = Neanderthal · light blue = Denisovan · taller = on both chromosomes · right = segment count');
});

/* ── mutation spectrum ── */
reveal('spectrum',(s,c)=>{
  const a=AC(), SUB=['C>A','C>G','C>T','T>A','T>C','T>G'];
  const COLS=[c.data2,c.ink,a.s,c.faintdata,a.n,c.data];
  const x0=44,x1=880,y0=30,H=190, data=D.spectrum, max=Math.max(...data.map(d=>d.frac)), bw=(x1-x0)/96, gap=bw*0.22;
  txt(s,{x:x0,y:14,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'占全部单碱基替换的比例':'SHARE OF ALL SINGLE-BASE SUBSTITUTIONS');
  data.forEach((d,i)=>{const gi=SUB.indexOf(d.sub), x=x0+i*bw, h=H*d.frac/max;
    const r=el(s,'rect',{x:x+gap/2,y:y0+H-h,width:bw-gap,height:h,rx:1,fill:COLS[gi]});
    tip(r,`${d.ctx}: ${(d.frac*100).toFixed(2)}% (${fmt(d.n)})`)});
  el(s,'line',{x1:x0,y1:y0+H,x2:x1,y2:y0+H,stroke:c.faint,'stroke-width':.9});
  SUB.forEach((sb,gi)=>{const xa=x0+gi*16*bw, xb=xa+16*bw;
    el(s,'rect',{x:xa+1,y:y0-12,width:16*bw-2,height:7,rx:1.5,fill:COLS[gi]});
    txt(s,{x:(xa+xb)/2,y:y0+H+16,'text-anchor':'middle','font-size':11,'font-weight':800,fill:c.ink},sb);
    const tot=data.filter(d=>d.sub===sb).reduce((p,q)=>p+q.frac,0);
    txt(s,{x:(xa+xb)/2,y:y0+H+30,'text-anchor':'middle','font-size':9,fill:c.muted},(tot*100).toFixed(1)+'%')});
  const ci=data.map((d,i)=>/\[C>T\]G/.test(d.ctx)?i:-1).filter(i=>i>=0);
  ci.forEach(i=>el(s,'circle',{cx:x0+i*bw+bw/2,cy:y0+H+42,r:2.2,fill:a.s}));
  txt(s,{x:x0+ci[0]*bw-6,y:y0+H+46,'text-anchor':'end','font-size':9,fill:c.muted},zh()?`CpG 位点的 C>T：${F.cpg}%`:`C>T at CpG: ${F.cpg}%`);
  foot(s,c,900,296,zh()?'96 根柱 = 6 种替换 × 前后各一个碱基的 16 种组合':'96 bars = 6 substitution types × 16 flanking-base combinations');
});

/* ── KIR ── */
reveal('kir',(s,c)=>{
  const a=AC();
  const ORDER=['KIR3DL3','KIR2DS2','KIR2DL2','KIR2DL3','KIR2DP1','KIR2DL1','KIR3DP1','KIR2DL4','KIR3DL1','KIR3DS1','KIR2DL5A','KIR2DL5B','KIR2DS3','KIR2DS5','KIR2DS1','KIR2DS4','KIR3DL2'];
  const ACT=new Set(['KIR2DS1','KIR2DS2','KIR2DS3','KIR2DS4','KIR2DS5','KIR3DS1']);
  const map=Object.fromEntries(D.kir.map(k=>[k.gene,k]));
  const x0=30,x1=880,y=64,bw=(x1-x0)/ORDER.length;
  txt(s,{x:x0,y:18,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'KIR 基因区（19q13.4）按染色体顺序':'THE KIR LOCUS (19q13.4) IN CHROMOSOMAL ORDER');
  el(s,'line',{x1:x0,y1:y+22,x2:x1,y2:y+22,stroke:c.quiet,'stroke-width':1});
  ORDER.forEach((g,i)=>{const k=map[g], on=k&&k.present, x=x0+i*bw;
    const r=el(s,'rect',{x:x+2,y:y,width:bw-4,height:44,rx:2.5,fill:on?(ACT.has(g)?a.s:a.n):'none',opacity:on?.9:1,
      stroke:on?'none':c.faint,'stroke-width':1.2,'stroke-dasharray':on?'':'3 3'});
    tip(r,on?`${g}: ${k.call}`:`${g}: absent`);
    txt(s,{x:x+bw/2,y:y-8,'text-anchor':'middle','font-size':8.5,'font-weight':on?700:500,fill:on?c.ink:c.faint,transform:`rotate(-42 ${x+bw/2} ${y-8})`},g.replace('KIR',''));
    if(on) txt(s,{x:x+bw/2,y:y+27,'text-anchor':'middle','font-size':9,fill:c.bg,'font-weight':800},'✓')});
  [[a.n,zh()?'抑制性受体':'inhibitory'],[a.s,zh()?'活化性受体':'activating'],['none',zh()?'缺失':'absent']].forEach(([col,lab],i)=>{
    const lx=x0+i*160;
    if(col==='none') el(s,'rect',{x:lx,y:y+44,width:11,height:11,rx:2,fill:'none',stroke:c.faint,'stroke-width':1.2,'stroke-dasharray':'3 3'});
    else el(s,'rect',{x:lx,y:y+44,width:11,height:11,rx:2,fill:col});
    txt(s,{x:lx+16,y:y+53,'font-size':9.5,fill:c.lab},lab)});
  txt(s,{x:x0,y:y+82,'font-size':12,'font-weight':800,fill:c.ink},zh()?'单倍型：AA（两条 A 单倍型）':'Haplotype: AA (two A haplotypes)');
  txt(s,{x:x0,y:y+99,'font-size':10,fill:c.lab},zh()?'HLA 配体：C*03:04 与 C*01:02 均为 C1 组 · B*13:01 为 Bw4 · B*54:01 为 Bw6':'HLA ligands: C*03:04 and C*01:02 are both C1 · B*13:01 is Bw4 · B*54:01 is Bw6');
});

/* ── stat grids and tables ── */
function statGrid(id,rows){const n=document.getElementById(id); if(!n)return;
  n.innerHTML=rows.map(([v,r,l])=>`<div class="kpi"><div class="v">${v}</div><div class="r">${r}</div><div class="l">${l}</div></div>`).join('')}
function deepRender(){
  const hd=document.getElementById('tb_hlad');
  if(hd) hd.innerHTML=D.hla_disease.map(r=>{const yes=r.carried==='yes';
    const cond=zh()?r.condition:r.condition_en, req=zh()?r.requirement:r.requirement_en, nt=zh()?r.note:r.note_en;
    return `<tr><td class="gene">${cond}</td><td class="rs">${req}</td><td class="g" style="color:${yes?'var(--ancs)':'var(--muted)'}">${yes?(zh()?'携带':'yes'):(zh()?'未携带':'no')}</td><td>${nt}</td></tr>`}).join('');
  const cd=document.getElementById('tb_cand');
  if(cd) cd.innerHTML=D.candidate.map(r=>`<tr><td class="gene">${r.gene}</td><td class="rs">${r.rsid} · ${zh()?r.variant:(r.variant_en||r.variant)}</td><td class="g">${r.genotype}</td><td class="rs">${zh()?r.popular_claim:(r.popular_claim_en||'')}</td><td>${zh()?r.what_evidence_supports:(r.evidence_en||'')}</td></tr>`).join('');
  statGrid('k_la',[[`${F.north}<small>%</small>`,zh()?`北方东亚成分，比北京汉平均高 ${F.sdchb} 个标准差`:`northern East Asian; ${F.sdchb} sd above the Beijing Han mean`,'northern east asian'],
                   [`${F.south}<small>%</small>`,zh()?`南方东亚成分（北京汉 ${100-F.chb}%、南方汉 ${100-F.chs}%）`:`southern (Beijing Han ${100-F.chb}%, southern Han ${100-F.chs}%)`,'southern east asian'],
                   [`${F.noise}<small>%</small>`,zh()?'欧洲 + 南亚，方法的噪声底':'European + South Asian, the noise floor','noise floor'],
                   [`${D.la_segments.filter(g=>g.anc==='SouthEA').length}`,zh()?'≥0.5 Mb 的南方片段，最长 6.6 Mb':'southern segments ≥0.5 Mb, longest 6.6 Mb','segments']]);
  statGrid('k_arch',[[`${F.archmb}<small>Mb</small>`,zh()?`古老人类片段总长，占常染色体 ${F.archpct}%`:`archaic span, ${F.archpct}% of the autosomes`,'archaic span'],
                     [`${F.archn}`,zh()?'片段数，其中 13 段为纯合':'segments, 13 of them homozygous','segments'],
                     [`${F.neamb}<small>Mb</small>`,zh()?'尼安德特来源':'from Neanderthals','neanderthal'],
                     [`${F.denmb}<small>Mb</small>`,zh()?'丹尼索瓦来源':'from Denisovans','denisovan']]);
  statGrid('k_phase',[[`${F.phpct}<small>%</small>`,zh()?`杂合位点已定相（${F.phased} / ${F.phhet}）`:`heterozygous sites phased (${F.phased} of ${F.phhet})`,'phased'],
                      [`${F.n50}<small>Mb</small>`,zh()?`相位块 N50，共 ${F.blocks} 个块`:`phase block N50, ${F.blocks} blocks`,'block n50'],
                      [`2`,zh()?'因此确定的顺反关系：UGT1A1 与 FUT3':'cis/trans questions settled: UGT1A1 and FUT3','cis / trans']]);
  statGrid('k_som',[[`${F.mtcn}`,zh()?'线粒体拷贝/细胞（正常 100–500）':'mtDNA copies per cell (100–500 normal)','mtdna copy number'],
                    [`${F.ydr}`,zh()?'Y 染色体深度比（<0.45 才提示镶嵌性丢失）':'Y depth ratio (mosaic loss below 0.45)','y depth ratio'],
                    [`0`,zh()?`19 个克隆性造血热点的支持读段，中位深度 ${D.chip_hotspots.median_depth}×`:`reads at 19 clonal-haematopoiesis hotspots, median depth ${D.chip_hotspots.median_depth}×`,'clonal haematopoiesis'],
                    [`${F.telk7}`,zh()?'端粒重复读段（6.52 亿条中），不足以估计端粒长度':'telomeric reads of 652 million, too few for a length estimate','telomeric reads']]);
  const tb=document.getElementById('tb_anc');
  if(tb) tb.innerHTML=D.ho_near_ancient.slice(0,10).map(r=>{const nm=zh()?(ANC_ZH[r.label]||r.label):r.label.replace(/_/g,' ');
    return `<tr><td class="gene">${nm}</td><td class="num">${r.n}</td><td class="num">${fmt(Math.round(r.date))}</td><td class="num">${r.d.toFixed(4)}</td></tr>`}).join('');
  const bb=document.getElementById('tb_blood');
  if(bb) bb.innerHTML=BLOODV3.map(r=>`<tr><td class="gene">${zh()?r.sys_zh:r.sys_en}</td><td class="g">${zh()?r.geno_zh:r.geno_en}</td><td class="gene">${zh()?r.ph_zh:r.ph_en}</td><td>${zh()?r.ev_zh:r.ev_en}</td></tr>`).join('');
}

/* ── circular genome overview ── */
reveal('circos',(s,c)=>{
  const a=AC(), CX=430, CY=430, GAP=0.9*Math.PI/180;
  const total=CHR22.reduce((p,ch)=>p+D.chrlen[ch],0);
  const span=2*Math.PI-CHR22.length*GAP;
  let ang=-Math.PI/2+GAP/2; const arcs={};
  CHR22.forEach(ch=>{const a0=ang,a1=ang+span*D.chrlen[ch]/total; arcs[ch]=[a0,a1]; ang=a1+GAP});
  const P=(r,t)=>[CX+r*Math.cos(t),CY+r*Math.sin(t)];
  const arc=(r0,r1,a0,a1,fill,op)=>{
    const [x0,y0]=P(r1,a0),[x1,y1]=P(r1,a1),[x2,y2]=P(r0,a1),[x3,y3]=P(r0,a0);
    const big=(a1-a0)>Math.PI?1:0;
    return el(s,'path',{d:`M${x0} ${y0}A${r1} ${r1} 0 ${big} 1 ${x1} ${y1}L${x2} ${y2}A${r0} ${r0} 0 ${big} 0 ${x3} ${y3}Z`,fill,opacity:op===undefined?1:op});
  };
  const at=(ch,pos)=>{const [a0,a1]=arcs[ch];return a0+(a1-a0)*pos/D.chrlen[ch]};
  const R={lab:410, ideo:[372,392], dens:[318,368], la:[292,314], arch:[272,288], roh:[256,268]};
  CHR22.forEach(ch=>{
    const [a0,a1]=arcs[ch];
    arc(R.ideo[0],R.ideo[1],a0,a1,c.track);
    const cen=D.cen[ch];
    if(cen){const t=at(ch,cen*1e6); const [x0,y0]=P(R.ideo[0],t),[x1,y1]=P(R.ideo[1],t);
      el(s,'line',{x1:x0,y1:y0,x2:x1,y2:y1,stroke:c.faint,'stroke-width':1.4})}
    const am=(a0+a1)/2, [lx,ly]=P(R.lab,am);
    txt(s,{x:lx,y:ly,'text-anchor':'middle','dominant-baseline':'central','font-size':11,'font-weight':700,fill:c.lab},ch);
    for(let p=0;p<D.chrlen[ch];p+=50e6){const t=at(ch,p);const [x0,y0]=P(R.ideo[1],t),[x1,y1]=P(R.ideo[1]+4,t);
      el(s,'line',{x1:x0,y1:y0,x2:x1,y2:y1,stroke:c.faint,'stroke-width':.7})}
  });
  const dens=Object.fromEntries(D.density.map(d=>[d.chrom,d.counts]));
  const allc=D.density.flatMap(d=>d.counts.filter(x=>x>0)).sort((p,q)=>p-q);
  const lo=allc[Math.floor(allc.length*0.02)], hi=allc[Math.floor(allc.length*0.98)];
  el(s,'circle',{cx:CX,cy:CY,r:R.dens[0],fill:'none',stroke:c.grid,'stroke-width':.7});
  el(s,'circle',{cx:CX,cy:CY,r:R.dens[1],fill:'none',stroke:c.grid,'stroke-width':.7});
  CHR22.forEach(ch=>{
    const co=dens[ch]||[]; if(co.length<2) return;
    let d=`M${P(R.dens[0],at(ch,0))[0]} ${P(R.dens[0],at(ch,0))[1]}`;
    co.forEach((v,i)=>{const t=at(ch,Math.min(i*5e6+2.5e6,D.chrlen[ch]));
      const f=Math.max(0,Math.min(1,(v-lo)/(hi-lo)));
      const p=P(R.dens[0]+(R.dens[1]-R.dens[0])*f,t); d+=`L${p[0].toFixed(1)} ${p[1].toFixed(1)}`});
    const [ex,ey]=P(R.dens[0],at(ch,D.chrlen[ch])); d+=`L${ex} ${ey}Z`;
    el(s,'path',{d,fill:c.data2,opacity:.34});
  });
  CHR22.forEach(ch=>{const [a0,a1]=arcs[ch]; arc(R.la[0],R.la[1],a0,a1,a.n,.5)});
  D.la_segments.forEach(g=>{
    const col=g.anc==='SouthEA'?a.s:c.fd;
    const a0=at(g.chrom,g.start),a1=at(g.chrom,g.end);
    const h=(R.la[1]-R.la[0])/2, r0=g.hap===1?R.la[0]+h:R.la[0];
    const p=arc(r0,r0+h,a0,Math.max(a1,a0+0.0016),col,.95);
    tip(p,`chr${g.chrom}:${fmt(g.start)}-${fmt(g.end)} · ${g.anc} · hap ${g.hap}`);
  });
  CHR22.forEach(ch=>{const [a0,a1]=arcs[ch]; arc(R.arch[0],R.arch[1],a0,a1,c.track)});
  D.archaic.forEach(g=>{
    const col=g.src==='Denisovan'?c.data2:a.a;
    const p=arc(R.arch[0],R.arch[1],at(g.chrom,g.start),Math.max(at(g.chrom,g.end),at(g.chrom,g.start)+0.0016),col,g.hom?1:.7);
    tip(p,`chr${g.chrom}:${fmt(g.start)}-${fmt(g.end)} · ${g.src} · ${g.mb} Mb${g.hom?' · homozygous':''}`);
  });
  CHR22.forEach(ch=>{const [a0,a1]=arcs[ch]; arc(R.roh[0],R.roh[1],a0,a1,c.track)});
  D.roh.forEach(r=>{ if(!arcs[r.chrom]) return;
    const p=arc(R.roh[0],R.roh[1],at(r.chrom,r.start),Math.max(at(r.chrom,r.end),at(r.chrom,r.start)+0.0016),c.ink,r.mb>=1.5?.95:.45);
    tip(p,`chr${r.chrom}:${fmt(r.start)}-${fmt(r.end)} · ${r.mb} Mb`);
  });
  const DEL=[["1",110230156,"GSTM1",14],["1",152555540,"LCE3B/C",34],["19",41349539,"CYP2A6",34],["19",54800853,"LILRA3",14],["22",24274142,"GSTT1",54],["4",69373820,"UGT2B17",14]];
  DEL.forEach(([ch,pos,name,drop])=>{
    const t=at(ch,pos); const [x0,y0]=P(R.roh[0]-3,t),[x1,y1]=P(R.roh[0]-drop,t);
    el(s,'line',{x1:x0,y1:y0,x2:x1,y2:y1,stroke:c.faint,'stroke-width':.8});
    const right=Math.cos(t)>=0;
    txt(s,{x:x1+(right?4:-4),y:y1,'text-anchor':right?'start':'end','dominant-baseline':'central','font-size':9,'font-weight':700,fill:c.lab,stroke:c.bg,'stroke-width':2.6,'paint-order':'stroke'},name);
  });
  txt(s,{x:CX,y:CY-16,'text-anchor':'middle','font-size':12,'font-weight':700,fill:c.muted,'letter-spacing':'.16em'},NAME().toUpperCase());
  txt(s,{x:CX,y:CY+14,'text-anchor':'middle','font-size':30,'font-weight':800,fill:c.ink,...NUM},'2.87 Gb');
  txt(s,{x:CX,y:CY+34,'text-anchor':'middle','font-size':10,fill:c.muted},zh()?'22 条常染色体':'22 autosomes');
  const RL=[[R.dens,zh()?'变异密度':'variant density',c.data2,.34],
            [R.la,zh()?'局部祖源':'local ancestry',a.n,.5],
            [R.arch,zh()?'古老人类片段':'archaic segments',a.a,1],
            [R.roh,zh()?'纯合区':'homozygous runs',c.ink,.95]];
  const kx=6, ky=18;
  RL.forEach(([rr,lab,col,op],i)=>{const y=ky+i*18;
    el(s,'rect',{x:kx,y:y-6,width:24,height:9,rx:1.5,fill:c.track});
    el(s,'rect',{x:kx+(i===0?0:5),y:y-6,width:i===0?24:11,height:9,rx:1.5,fill:col,opacity:op});
    txt(s,{x:kx+30,y:y+2,'font-size':9.5,fill:c.lab},lab)});
  const yo=ky+RL.length*18;
  el(s,'rect',{x:kx,y:yo-6,width:24,height:9,rx:1.5,fill:a.n,opacity:.5});
  el(s,'rect',{x:kx+5,y:yo-6,width:11,height:9,rx:1.5,fill:a.s});
  txt(s,{x:kx+30,y:yo+2,'font-size':9.5,fill:c.lab},zh()?'南方东亚片段':'southern segment');
});

/* ── archaic gene families ── */
reveal('afam',(s,c)=>{
  const a=AC(), rows=D.archaic_families.slice(0,10), x0=240,x1=760,top=16,rh=25;
  const X=v=>x0+(x1-x0)*v/100;
  [0,25,50,75,100].forEach(v=>{el(s,'line',{x1:X(v),y1:top-8,x2:X(v),y2:top+rows.length*rh-6,stroke:c.grid,'stroke-width':.6});
    txt(s,{x:X(v),y:top+rows.length*rh+8,'text-anchor':'middle','font-size':8,fill:c.faint},v+'%')});
  txt(s,{x:x0,y:8,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'该基因家族有多少比例的基因落在渗入片段内':'SHARE OF THE FAMILY\'S GENES THAT FALL INSIDE AN INTROGRESSED SEGMENT');
  rows.forEach((r,i)=>{const y=top+i*rh+6, w=X(r.pct)-x0, hi=r.pct>=25;
    txt(s,{x:x0-10,y:y+3,'text-anchor':'end','font-size':10,'font-weight':hi?800:600,fill:hi?c.ink:c.lab},zh()?r.family:(r.family_en||r.family));
    const b=el(s,'rect',{x:x0,y:y-6,width:Math.max(w,1.5),height:12,rx:2,fill:hi?a.a:c.data2,opacity:hi?1:.7,class:'fade',style:`animation-delay:${i*.05}s`});
    tip(b,`${zh()?r.family:(r.family_en||r.family)}: ${r.carried}/${r.total} · ${r.examples}`);
    txt(s,{x:X(r.pct)+7,y:y+3,'font-size':9.5,'font-weight':800,fill:c.ink},`${r.carried}/${r.total}`)});
  foot(s,c,900,296,zh()?`片段总共覆盖 ${D.archaic_genes.n_genes} 个蛋白编码基因，随机打乱片段位置的期望是 ${Math.round(D.archaic_genes.null_mean)} 个：总数不高于随机，特定家族才是重点`
                     :`the segments cover ${D.archaic_genes.n_genes} protein-coding genes against ${Math.round(D.archaic_genes.null_mean)} expected by chance: the total is unremarkable, the families are the point`);
});

/* ── behaviour PRS ── */
reveal('behaviour',(s,c)=>{
  const a=AC(), rows=[...D.behaviour].sort((p,q)=>q.pct_EAS-p.pct_EAS), x0=230,x1=800,top=24,rh=23;
  const X=v=>x0+(x1-x0)*v/100;
  el(s,'rect',{x:X(25),y:top-10,width:X(75)-X(25),height:rows.length*rh+2,fill:c.track,opacity:.6});
  [0,25,50,75,100].forEach(v=>{el(s,'line',{x1:X(v),y1:top-10,x2:X(v),y2:top+rows.length*rh-4,stroke:v===50?c.faint:c.grid,'stroke-width':v===50?.9:.5});
    txt(s,{x:X(v),y:top+rows.length*rh+10,'text-anchor':'middle','font-size':8,fill:c.faint},v)});
  txt(s,{x:x0,y:10,'font-size':8,'font-weight':600,fill:c.faint,'letter-spacing':'.08em'},zh()?'在 504 名东亚参考个体中的百分位':'PERCENTILE AMONG THE 504 EAST ASIAN REFERENCE INDIVIDUALS');
  rows.forEach((r,i)=>{const y=top+i*rh+6, eas=r.panel==='EAS';
    txt(s,{x:x0-12,y:y+3,'text-anchor':'end','font-size':10,'font-weight':eas?700:500,fill:eas?c.ink:c.muted},zh()?r.zh:r.en);
    const mx=X(r.pct_EAS);
    if(eas){const d=el(s,'rect',{x:mx-4.5,y:y-4.5,width:9,height:9,fill:a.n,stroke:c.bg,'stroke-width':1.2,transform:`rotate(45 ${mx} ${y})`,class:'pop'});
      d.style.animationDelay=(i*25)+'ms'; tip(d,`${r.en}: ${r.pct_EAS} pct · z ${r.z} · coverage ${r.coverage}%`)}
    else {const d=el(s,'circle',{cx:mx,cy:y,r:4.5,fill:c.bg,stroke:a.s,'stroke-width':1.8,class:'pop'});
      d.style.animationDelay=(i*25)+'ms'; tip(d,`${r.en}: ${r.pct_EAS} pct · z ${r.z} · coverage ${r.coverage}% · European-derived`)}
    txt(s,{x:x1+12,y:y+3,'font-size':9.5,'font-weight':800,fill:eas?c.ink:c.muted},r.pct_EAS.toFixed(0));
    txt(s,{x:x1+36,y:y+3,'font-size':8,fill:c.faint},r.coverage+'%')});
  txt(s,{x:x1+12,y:top-10,'font-size':7.5,'font-weight':600,fill:c.faint,'letter-spacing':'.06em'},zh()?'百分位 · 覆盖':'PCT · COV');
  foot(s,c,900,426,zh()?'菱形 = 东亚人群训练 · 空心圆 = 欧洲人群训练 · 灰带 = 中间一半的人':'diamond = trained in East Asians · hollow circle = trained in Europeans · band = middle half');
});

renderAll();
