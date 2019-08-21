let SessionLoad = 1
if &cp | set nocp | endif
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
silent tabonly
cd /project/gilad/bjf79/projects/Comparative_eQTL/code/snakemake_workflow
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
argglobal
%argdel
$argadd Snakefile
edit rules/eQTL_analysis.smk
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
argglobal
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let s:l = 318 - ((27 * winheight(0) + 25) / 50)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
318
normal! 09|
tabnext 1
badd +8 Snakefile
badd +313 rules/eQTL_analysis.smk
badd +76 rules/RNASeqMapping.smk
badd +259 rules/populationstructure.smk
badd +527 rules/eqtl_calling.smk
badd +49 /project/gilad/bjf79/software/locusZoomWithLD/locuszoom/conf/m2zfast.conf
badd +112 cluster-config.json
badd +101 scripts/MatrixEqtl_Cis.AllPvals.R
badd +128 rules/PowerAnalysis.smk
badd +154 config.yaml
badd +66 rules/sqtl_calling.smk
badd +7 scripts/MakeMetaPlot.sh
badd +1 .gitignore
badd +24 scratch/TODO.txt
badd +7 scripts/StandardizeAndQuantileNormalize.py
badd +79 rules/calling.smk
badd +8 rules/common.smk
badd +38 ~/CurrentProjects/recurrent-splicing-identification/rules/leafcutter.smk
badd +1 ~/CurrentProjects/201905_MPESeq/rules/STAR_Alignment.smk
badd +12 snakemake.sbatch
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToO
set winminheight=1 winminwidth=1
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
nohlsearch
let g:this_session = v:this_session
let g:this_obsession = v:this_session
let g:this_obsession_status = 2
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
