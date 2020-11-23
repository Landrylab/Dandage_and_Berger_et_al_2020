from rohan.global_imports import *
import networkx as nx
def plot_network(dplot,params,ax=None,title=None,test=False):
    """
    Plot the PPI networks within complexes
    """
    g=nx.from_pandas_edgelist(dplot, **params['params_from_pandas_edgelist'])
    pos=nx.circular_layout(g)
    for nodetype in params['nodetype2draw_networkx_nodes']:
        nodes=nx.draw_networkx_nodes(g,pos,
                               **params['nodetype2draw_networkx_nodes'][nodetype],
                               ax=ax,
                              )
        if not nodes is None:
            nodes.set_edgecolor(params['nodetype2draw_networkx_nodes'][nodetype]['node_color'])
            nodes.set_facecolor('w')

    edge2width=nx.get_edge_attributes(g,'interaction score')
    color2edgelist=groupby_value(nx.get_edge_attributes(g,'edge color'))
    for color in color2edgelist:
        edgelist=color2edgelist[color]
        nx.draw_networkx_edges(g,pos,
                               edgelist=edgelist,
                               edge_color=color,
                               alpha=0.5  if color!='#428774' else 0.5,
                               style='solid' if color!='#428774' else 'dashed',
                               width=[edge2width[e]*5 for e in edgelist],
                                     ax=ax,
                              )
    from rohan.dandage.plot.colors import saturate_color,rgb2hex,rgbfloat2int
    for nodetype in params['nodetype2draw_networkx_nodes']:
        nx.draw_networkx_labels(g,pos,
                    labels={k:k.split(' ')[-1] for k in params['nodetype2draw_networkx_nodes'][nodetype]['nodelist']},
                    font_color=rgb2hex(rgbfloat2int(saturate_color(params['nodetype2draw_networkx_nodes'][nodetype]['node_color'],1.2))),
                                font_weight='bold',
                                ax=ax,
                               )        
    if hasattr(dplot,'name'):
        ax.set_title(dplot.name)
    elif not title is None:
        ax.set_title(title)        
    ax.set_xlim(-1.25,1.25)
    ax.set_ylim(-1.25,1.25)
    if not test:
        ax.set_axis_off()
    else:
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    return ax
def sort_nodes_in_dplot(dplot):
    """
    Sort the nodes in the network by their names.
    """
    dplot.loc[:,'sort by']=dplot['interaction id'].apply(lambda x: x.split(' ')[1])
    dplot=dplot.sort_values(by=['sort by'])
    df_=pd.DataFrame({'species name interactor':unique(dplot['species name interactor1'].tolist()+dplot['species name interactor2'].tolist())})
    df_=df_.sort_values('species name interactor')
    df_.loc[:,'gene name']=df_['species name interactor'].apply(lambda x: x.split(' ')[-1])
    genename2interactors=df_.groupby('gene name').agg(list).to_dict()['species name interactor']
    genename_sorted=list(sort_dict({k:[len(genename2interactors[k])] for k in genename2interactors},1).keys())[::-1]
    genename2interactors={k:genename2interactors[k] for k in [k_ for k_ in genename2interactors if len(genename2interactors[k_])==2]}
    df_=({'':[sorted(genename2interactors[k]) for k in genename2interactors]})
    df_=pd.DataFrame([sorted(genename2interactors[k]) for k in genename2interactors])
    if len(df_)!=0:
        df_.columns=[f'species name interactor{i}' for i in [1,2]]
        dplot=df_.append(dplot,sort=True)
    return dplot

def plot_interactome_hybrid(dplot,title=None,
                            ax=None,
                            plot=True,
                            test=False):
    """
    Plot PPI networks within the protein complexes of hybrid.
    """
    if not 'species name interactor1' in dplot:
        dplot=dplot.join(dplot['species name interaction id'].str.split('--').apply(pd.Series).rename(columns={i:f"species name interactor{int(i+1)}" for i in [0,1]}))
    if not 'species1 name' in dplot:
        dplot=dplot.join(dplot['species name interaction id'].apply(lambda x: [s.split(' ')[0] for s in x.split('--')]).apply(pd.Series).rename(columns={i:f"species{i+1} name" for i in [0,1]}))

    dplot=sort_nodes_in_dplot(dplot)
    l=[dplot.set_index(f'species name interactor{i}')[f'protein abundance (log10 scale) interactor{i}'].fillna(0).to_dict() for i in [1,2]]
    d=merge_dict(*l)
    node2size={k:d[k][0] for k in d}
    if test:
        pprint(node2size)
        l=np.random.randn(len(node2size.keys()))
        l=sorted((l-l.min())/(l.max()-l.min())*2)
        node2size={k:l[i] for i,k in enumerate(node2size.keys())}
        pprint(node2size)        
    node2size={k:scale_node_size(node2size[k]) for k in node2size}
    
    dplot['edge color']=dplot.apply(lambda x: np.nan if pd.isnull(x['species1 name']) else '#6cacc5' if x['species1 name']==x['species2 name'] and x['species1 name']=='Suva' else '#dc9f4e' if x['species1 name']==x['species2 name'] and x['species1 name']=='Scer' else '#428774' if x['species1 name']!=x['species2 name'] else np.nan,axis=1)
    nodetype2params={'Scer':{'node_color':'#dc9f4e'},
                   'Suva':{'node_color':'#6cacc5'},
                    }
    for nodetype in nodetype2params:
        nodetype2params[nodetype]['nodelist']=[s for s in node2size if s.startswith(nodetype)] 
        nodetype2params[nodetype]['node_size']=[node2size[n] for n in nodetype2params[nodetype]['nodelist']]
        nodetype2params[nodetype]['zorder']=-3
    
    params={'params_from_pandas_edgelist':{'source':'species name interactor1',
                                           'target':'species name interactor2',
                                         'edge_attr':[
                                                     'interaction score', #colorby
                                                     'edge color',
                                                     ], 
                                          },
            'draw_networkx_edges':{'edge_vmin':0,'edge_vmax':1,'edge_color':'lightgray'},
            'nodetype2draw_networkx_nodes':nodetype2params,
           }
    if test:
        to_table(dplot,"test/dplot.tsv")
        to_dict(params,"test/params.yml")
    if plot:
        return plot_network(dplot,params,test=test,ax=ax)
    else:
        g=nx.from_pandas_edgelist(dplot, **params['params_from_pandas_edgelist'])
        return g
    
def plot_complex_hybrid_parent(funs):
    """
    Plot PPI networks within the protein complexes of parental species.
    """
    complexid2size=read_table(f"{dirname(__file__)}/plot/circle_complexes_hybrid.tsv").set_index('complex id')['complex size'].to_dict()
    complexid=funs[0].__name__.split('_')[-1]
    plt.figure(figsize=[8,(complexid2size[f"CPX-{complexid}"]*1.7)+2])
    axs=[plt.subplot2grid([2,1],[0,0],1,1),
         plt.subplot2grid([2,1],[1,0],1,1)
        ]
    for fun,ax in zip(funs,axs):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            ax=fun(ax=ax)
            ax.set_title(fun.__name__.split('_')[-3],loc='left')
    labelsubplots(axes=np.ravel(axs),xoff=0)
    plt.tight_layout()
from rohan.global_imports import *

## FigureS01:panel0
def plot_bar_peptides_detected_by_unique_not_unique(plotp="plot/bar_peptides_detected_by_unique_not_unique.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    dplot.set_index('subset\n(total peptides)').loc[:,['unique peptides','non-unique peptides']].plot.barh(stacked=True,ax=ax)
    ax.legend(bbox_to_anchor=[1,1])
    ax.set(xlim=[0,100],xlabel='%')
    dplot['y']=range(len(dplot))
    dplot.apply(lambda x: ax.text(1,x['y'],f"{x['unique peptides']:.0f}%",color='darkred'),axis=1)
    dplot.apply(lambda x: ax.text(100,x['y'],f"{x['non-unique peptides']:.0f}%",color='gray',ha='right'),axis=1)
    return ax


## FigureS01:panel1
def plot_bar_proteins_detected_by_unique_not_unique(plotp="plot/bar_proteins_detected_by_unique_not_unique.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    dplot.set_index('subset\n(total proteins detected)').loc[:,['detected by\nunique peptides','detected\nexclusively by\nnon-unique peptides']].plot.barh(stacked=True,ax=ax)
    ax.legend(bbox_to_anchor=[1,1])
    ax.set(xlim=[0,100],xlabel='%')
    dplot['y']=range(len(dplot))
    dplot.apply(lambda x: ax.text(1,x['y'],"{:.0f}%".format(x['detected by\nunique peptides']),
    color='darkred'),axis=1)
    dplot.apply(lambda x: ax.text(100,x['y'],"{:.0f}%".format(x['detected\nexclusively by\nnon-unique peptides']),
    color='gray',ha='right'),axis=1)
    return ax


## FigureS01:panel2
def plot_hist_orthologs_dissimilar_peptides_per_protein__Jaccard_distance___(plotp="plot/hist_orthologs_dissimilar_peptides_per_protein__Jaccard_distance___.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=sns.violinplot(data=dplot,y='proteins',x=params['x'],palette=['gray','#d24043'],
    scale='width',ax=ax)
    [ax.text(params['label2x'][s.get_text()],si-0.1,params['label2x'][s.get_text()],ha='center',color='k') for si,s in enumerate(ax.get_yticklabels())]
    ax.set_xlim(0,100)
    ax.set_ylabel('')
    return ax


## FigureS01:panel3
def plot_hist_orthologs_protein_sequence_divergence____(plotp="plot/hist_orthologs_protein_sequence_divergence____.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=sns.violinplot(data=dplot,y='proteins',x=params['x'],palette=['gray','#d24043'],
    scale='width',ax=ax)
    [ax.text(params['label2x'][s.get_text()],si-0.1,params['label2x'][s.get_text()],ha='center',color='k') for si,s in enumerate(ax.get_yticklabels())]
    ax.set_xlim(0,100)
    ax.set_ylabel('')
    return ax


## FigureS03:panel0
def plot_heatmap_interaction_qc_interaction_bool_db_ROC_AUC(plotp="plot/heatmap_interaction_qc_interaction_bool_db_ROC_AUC.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    from rohan.dandage.plot.annot import annot_heatmap
    if ax is None:ax=plt.subplot()
    dplot_=dplot.pivot_table(index='quantile threshold',columns='window size',values=params['colvalue'])
    ax=sns.heatmap(dplot_,
    cmap='Reds',cbar_kws={'label':params['label']},ax=ax)
    dannot=dplot_.applymap(lambda x : 'max' if dplot_.max().max()==x else '')
    dannot.columns=[str(c) for c in dannot.columns]
    dannot.index=[str(c) for c in dannot.index]
    annot_heatmap(ax, dannot.T, xoff=-0.1, yoff=0.2, kws_text={'color':'limegreen'}, annot_left='(', annot_right=')', annothalf='upper')
    from rohan.dandage.plot.ax_ import set_logo
    set_logo(imp=f'{dirname(__file__)}/var/logos/Scer.svg.png',size=0.2,ax=ax)
    return ax


## Figure02:panel0
def plot_line_interaction_qc_ROC_AUC_DTW__window_04__score__quantile_0_75__interaction_bool_db_formated(plotp="plot/line_interaction_qc_ROC_AUC_DTW__window_04__score__quantile_0_75__interaction_bool_db_formated.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    from rohan.dandage.plot.ax_ import set_equallim
    if ax is None:ax=plt.subplot()
    ax.plot(dplot['FPR'], dplot['TPR'], 
    color=params['color'],
    label=params['label'],lw=2)
    ax.plot([0, 1], [0, 1], color='gray',linestyle='--')
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xlabel('FPR')
    ax.set_ylabel('TPR')
    ax.legend()
    from rohan.dandage.plot.ax_ import format_ticklabels
    ax=format_ticklabels(ax)
    from rohan.dandage.plot.ax_ import set_logo
    set_logo(imp=f"{dirname(__file__)}/var/logos/Scer.svg",size=0.2,ax=ax)
    return ax


## Figure02:panel1
def plot_line_interaction_score_protein_dosage_balance_between_hybrid_Scer_Suva(plotp="plot/line_interaction_score_protein_dosage_balance_between_hybrid_Scer_Suva.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    dplot.plot(x=f"{params['colx']}\n(qbin mid-point)",y='interaction score mean',
    yerr='interaction score confidence_interval_95',color=params['color'],
    linestyle='--' if 'inter' in params['label'] else '-',
    lw=3,
    ax=ax,
    label=params['label'])
    ax.set(**{'ylabel':'interaction score'})
    from rohan.dandage.plot.ax_ import format_ticklabels
    ax=format_ticklabels(ax)
    return ax


## Figure02:panel2
def plot_line_interaction_score_protein_dosage_balance_within_hybrid_Scer_Scer(plotp="plot/line_interaction_score_protein_dosage_balance_within_hybrid_Scer_Scer.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    dplot.plot(x=f"{params['colx']}\n(qbin mid-point)",y='interaction score mean',
    yerr='interaction score confidence_interval_95',color=params['color'],
    linestyle='--' if 'inter' in params['label'] else '-',
    lw=3,
    ax=ax,
    label=params['label'])
    ax.set(**{'ylabel':'interaction score'})
    from rohan.dandage.plot.ax_ import format_ticklabels
    ax=format_ticklabels(ax)
    return ax


## Figure02:panel3
def plot_line_interaction_score_protein_dosage_balance_within_hybrid_Suva_Suva(plotp="plot/line_interaction_score_protein_dosage_balance_within_hybrid_Suva_Suva.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    dplot.plot(x=f"{params['colx']}\n(qbin mid-point)",y='interaction score mean',
    yerr='interaction score confidence_interval_95',color=params['color'],
    linestyle='--' if 'inter' in params['label'] else '-',
    lw=3,
    ax=ax,
    label=params['label'])
    ax.set(**{'ylabel':'interaction score'})
    from rohan.dandage.plot.ax_ import format_ticklabels
    ax=format_ticklabels(ax)
    return ax


## Figure02:panel4
def plot_line_interaction_score_protein_dosage_balance_within_parent_Scer_Scer(plotp="plot/line_interaction_score_protein_dosage_balance_within_parent_Scer_Scer.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    dplot.plot(x=f"{params['colx']}\n(qbin mid-point)",y='interaction score mean',
    yerr='interaction score confidence_interval_95',color=params['color'],
    linestyle='--' if 'inter' in params['label'] else '-',
    lw=3,
    ax=ax,
    label=params['label'])
    ax.set(**{'ylabel':'interaction score'})
    from rohan.dandage.plot.ax_ import format_ticklabels
    ax=format_ticklabels(ax)
    return ax


## Figure02:panel5
def plot_line_interaction_score_protein_dosage_balance_within_parent_Suva_Suva(plotp="plot/line_interaction_score_protein_dosage_balance_within_parent_Suva_Suva.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    dplot.plot(x=f"{params['colx']}\n(qbin mid-point)",y='interaction score mean',
    yerr='interaction score confidence_interval_95',color=params['color'],
    linestyle='--' if 'inter' in params['label'] else '-',
    lw=3,
    ax=ax,
    label=params['label'])
    ax.set(**{'ylabel':'interaction score'})
    from rohan.dandage.plot.ax_ import format_ticklabels
    ax=format_ticklabels(ax)
    return ax


## Figure02:panel6
def plot_scatter_interaction_score_distance_proteasome(plotp="plot/scatter_interaction_score_distance_proteasome.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from rohan.dandage.plot.scatter import plot_reg
    ax=plot_reg(dplot,'distance $\AA$','interaction score',
    scafmt='sca',
    title_stat=True,
    params_scatter={'color':'#d24043'},
    ax=ax,trendline_lowess=True)
    ax.scatter(dplot.loc[dplot['interaction bool biogrid direct interaction'],'distance $\AA$'], 
    dplot.loc[dplot['interaction bool biogrid direct interaction'],'interaction score'], 
    s=80, facecolors='none', edgecolors='k',label='direct PPIs')
    ax.legend(loc=1,frameon=True)
    ax.set_ylim(-0.05,1.05)
    from rohan.dandage.plot.ax_ import format_ticklabels
    ax=format_ticklabels(ax)
    from rohan.dandage.plot.ax_ import set_logo
    set_logo(imp=f"{dirname(__file__)}/var/logos/complex_proteasome_Scer.svg.png",size=0.2,ax=ax)
    return ax


## Figure02:panel7
def plot_scatter_protein_abundance_db(plotp="plot/scatter_protein_abundance_db.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from rohan.dandage.plot.scatter import plot_scatter
    ax=plot_scatter(dplot,
    colx='protein abundance\n(log10 scale)',
    coly='protein abundance\nppm (log10 scale)',
    colz=None,
    kind='hexbin',
    trendline_method='poly',
    stat_method="spearman",
    bootstrapped=False,
    params_plot={},
    cmap='Reds',label_colorbar=None,
    gridsize=25,
    params_plot_trendline={},
    params_set_label={'title':True},
    ax=ax,)
    ax.set_ylabel(ax.get_ylabel().replace('PAXdb ',''))
    from rohan.dandage.plot.ax_ import format_ticklabels
    ax=format_ticklabels(ax)
    from rohan.dandage.plot.ax_ import set_logo
    ax.set(xlim=[0.4,3.25],
    ylim=[1.25,3.75])
    set_logo(imp=f"{dirname(__file__)}/var/logos/Scer.svg.png",size=0.2,ax=ax)
    return ax


## Figure02:panel8
def plot_violin_interaction_score_mean_by_protein_complexes_and_random(plotp="plot/violin_interaction_score_mean_by_protein_complexes_and_random.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    from rohan.dandage.plot.dist import plot_dist_comparison
    params_ax={
    }        
    if ax is None:ax=plt.subplot()
    ax=plot_dist_comparison(dplot,ax=ax,**params,params_ax=params_ax)
    ax.set_ylabel('interaction score')
    ax.set_xlim(-0.45,0.45)
    from rohan.dandage.plot.ax_ import set_logo
    set_logo(imp=f'{dirname(__file__)}/var/logos/Scer.svg.png',size=0.2,ax=ax)
    return ax


## Figure02:panel9
def plot_violin_interaction_score_mean_by_protein_subcomplex(plotp="plot/violin_interaction_score_mean_by_protein_subcomplex.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    from rohan.dandage.plot.dist import plot_dist_comparison
    params_ax={
    'xlabel':'',
    }        
    if ax is None:ax=plt.subplot(111)
    ax=plot_dist_comparison(dplot,ax=ax,**params,params_ax=params_ax)
    ax.set_ylabel('interaction score')
    from rohan.dandage.plot.ax_ import set_logo
    set_logo(imp=f'{dirname(__file__)}/var/logos/complex_proteasome_Scer.svg',size=0.2,ax=ax)
    return ax


## FigureS04:panel0
def plot_network_within_parent_Scer_Scer_proteasome_CPX_2262(plotp="plot/network_within_parent_Scer_Scer_proteasome_CPX_2262.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    plot_network_legend(ax,legend_color=False)
    ax.set_axis_off()
    return ax


## FigureS05:panel0
def plot_line_interaction_score_protein_complexes_and_random(plotp="plot/line_interaction_score_protein_complexes_and_random.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=sns.pointplot(data=dplot,x=params['colx'],y=params['coly'],
    hue=params['colhue'],ci=95,dodge=True,
    order=params['xs'],hue_order=params['hues'],
    linestyles=params['linestyles'],
    palette=params['palette'],
    ax=ax)
    ax.legend(bbox_to_anchor=[1,1])
    ax.set_xlabel('')
    ax.set_xticklabels([t.get_text().replace(' ','\n') for t in ax.get_xticklabels()])
    return ax


## FigureS05:panel1
def plot_line_protein_dosage_balance_protein_complexes_and_random(plotp="plot/line_protein_dosage_balance_protein_complexes_and_random.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=sns.pointplot(data=dplot,x=params['colx'],y=params['coly'],
    hue=params['colhue'],ci=95,dodge=True,
    order=params['xs'],hue_order=params['hues'],
    linestyles=params['linestyles'],
    palette=params['palette'],
    ax=ax)
    ax.legend(bbox_to_anchor=[1,1])
    ax.set_xlabel('')
    ax.set_xticklabels([t.get_text().replace(' ','\n') for t in ax.get_xticklabels()])
    return ax


## Figure03:panel0
def plot_circle_complexes_hybrid(plotp="plot/circle_complexes_hybrid.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    from rohan.dandage.plot.scatter import plot_circlify
    from rohan.dandage.plot.colors import get_cmap_subset
    if ax is None:ax=plt.subplot()
    dplot=dplot.sort_values(by='median interaction score per protein',ascending=False)
    circvar2col={'parent id':'complex name',
    'parent datum':'complex size',
    'parent color':'idx',
    'child id':'species name interactor',
    }
    lineside2params=plot_circlify(dplot,circvar2col=circvar2col,threshold_side=-0.26,
    cmap_child=get_cmap_subset('binary', vmin=0.05, vmax=0.2, n=10),
    cmap_parent=None,#get_cmap_subset('binary', vmin=0.1, vmax=0.4, n=100),                              
    ax=ax,
    )
    return ax


## Figure03:panel1
def plot_contour_within_hybrids__Scer__Scer__within_hybrids__Suva__Suva__between_hybrids__Scer__Suva_(plotp="plot/contour_within_hybrids__Scer__Scer__within_hybrids__Suva__Suva__between_hybrids__Scer__Suva_.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    from rohan.dandage.plot.contour import plot_contourf
    from rohan.dandage.plot.colors import get_cmap_subset
    fig,ax=plot_contourf(
    x=dplot[params['colx']].values,
    y=dplot[params['coly']].values,
    z=dplot[params['colz']].values,
    grid_n=20,
    params_contourf={'cmap':get_cmap_subset('binary_r',0,0.8)},
    labelx=f"{params['colx']}",
    labely=f"{params['coly']}",
    labelz=f"{params['colz']}",
    params_cbar=params['params_cbar'],
    ax=ax,
    fig=fig,
    test=False)
    _=ax.text(*ax.get_ylim(),params['label'],va='bottom')
    return ax


## Figure03:panel2
def plot_contour_within_parents__Scer__Scer__within_parents__Suva__Suva__between_hybrids__Scer__Suva_(plotp="plot/contour_within_parents__Scer__Scer__within_parents__Suva__Suva__between_hybrids__Scer__Suva_.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    from rohan.dandage.plot.contour import plot_contourf
    from rohan.dandage.plot.colors import get_cmap_subset
    fig,ax=plot_contourf(
    x=dplot[params['colx']].values,
    y=dplot[params['coly']].values,
    z=dplot[params['colz']].values,
    grid_n=20,
    params_contourf={'cmap':get_cmap_subset('binary_r',0,0.8)},
    labelx=f"{params['colx']}",
    labely=f"{params['coly']}",
    labelz=f"{params['colz']}",
    params_cbar=params['params_cbar'],
    ax=ax,
    fig=fig,
    test=False)
    _=ax.text(*ax.get_ylim(),params['label'],va='bottom')
    return ax


## Figure03:panel3
def plot_line_connections_similarity_interaction_types_ms(plotp="plot/line_connections_similarity_interaction_types_ms.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from matplotlib.lines import Line2D
    legend_elements=[Line2D([0], [0], marker='o', color='none', label=f"{k} protein",markeredgewidth=4,
    markerfacecolor='none',markeredgecolor=params['element2color'][k],
    markersize=15) for k in ['Scer','Suva']]
    from rohan.dandage.plot.line import plot_connections
    ax=plot_connections(dplot=dplot,
    **{k:params[k] for k in params if k!='label2loc'},
    params_text={'ha':'center','va':'top'},
    ax=ax,test=False,legend_elements=legend_elements)
    label2loc=params['label2loc']
    _=[ax.add_patch(patch) for patch in [plt.Rectangle(**label2loc[k]) for k in label2loc]]
    _=[ax.text(label2loc[k]['xy'][0],
    label2loc[k]['xy'][1]+label2loc[k]['height']*0.5,k,color=label2loc[k]['edgecolor'],va='center',ha='center') for k in ['parents','hybrid']]
    _=[ax.text(label2loc[k]['xy'][0]+label2loc[k]['width']*0.5,
    label2loc[k]['xy'][1]+label2loc[k]['height'],k,color=label2loc[k]['edgecolor'],va='bottom',ha='center') for k in ['intralogous PPIs','interlogous PPIs']]
    return ax


## Figure03:panel4
def plot_scatter_interaction_score_interlogous_PPIs_predicted_ms(plotp="plot/scatter_interaction_score_interlogous_PPIs_predicted_ms.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from rohan.dandage.plot.ax_ import set_equallim
    from rohan.dandage.stat.corr import get_corr
    dplot.groupby('species type',sort=False).apply(lambda df: df.plot.scatter(x=params['colz'],
    y=f"{params['colz']} predicted",alpha=0.1,s=1,
    **{'color':params['element2color'][df.name],
    'label':'intralogous PPIs\nin '+df.name+'\n'+get_corr(df[params['colz']],
    df[f"{params['colz']} predicted"],
    method='spearman',
    bootstrapped=False,ci_type='max',outstr=True).replace('\n',' ')
    },
    ax=ax))
    _, labels = ax.get_legend_handles_labels()
    from matplotlib.lines import Line2D
    handles=[Line2D([0], [0], marker='o', color='none', 
    markerfacecolor=params['element2color'][k],
    markeredgecolor=params['element2color'][k],
    markersize=5) for k in ['parent' if 'parent' in s else 'hybrid' for s in labels]]
    _=ax.legend(handles=handles[::-1],labels=labels[::-1],
    title='predicted from',loc=2,bbox_to_anchor=[1,1])
    ax.axis('scaled')
    _=ax.set(**{'xlim':[0,1],'ylim':[0,1]})
    ax=set_equallim(ax,diagonal=True)
    return ax


## FigureS06:panel0
def plot_dendogram_interaction_score_within_species(plotp="plot/dendogram_interaction_score_within_species.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    dplot=dplot.set_index(['interaction id']) if 'interaction id' in dplot else dplot
    ax=sns.matrix.dendrogram(dplot,ax=ax)
    return ax


## FigureS06:panel1
def plot_dendogram_protein_abundance(plotp="plot/dendogram_protein_abundance.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    ax=sns.matrix.dendrogram(dplot.set_index(['interactor']),ax=ax)
    return ax


## FigureS06:panel2
def plot_heatmap_interaction_score_within_species(plotp="plot/heatmap_interaction_score_within_species.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    ax=sns.heatmap(dplot.set_index('interaction id'),cmap='Reds',
    cbar_kws={'label':'interaction score',"orientation": "horizontal",},
    ax=ax)
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_ylabel('PPIs')
    return ax


## FigureS06:panel3
def plot_heatmap_protein_abundance(plotp="plot/heatmap_protein_abundance.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    ax=sns.heatmap(dplot.set_index(['interactor']),cmap='Reds',
    cbar_kws={'label':'protein abundance\n(log10 scale)',"orientation": "horizontal",},
    ax=ax)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_ylabel('genes')
    return ax


## FigureS07:panel0
def plot_network_interactions_hybrid_CPX_2262(plotp="plot/network_interactions_hybrid_CPX_2262.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS07:panel1
def plot_network_interactions_parent_CPX_2262(plotp="plot/network_interactions_parent_CPX_2262.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS08:panel0
def plot_network_interactions_hybrid_CPX_1671(plotp="plot/network_interactions_hybrid_CPX_1671.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS08:panel1
def plot_network_interactions_parent_CPX_1671(plotp="plot/network_interactions_parent_CPX_1671.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS09:panel0
def plot_dist_interaction_type_interaction_score_complex_counts(plotp="plot/dist_interaction_type_interaction_score_complex_counts.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=sns.boxplot(data=dplot,**params['params_dist'],
    ax=ax,showbox=False,showcaps=False,showfliers=False)
    ax=sns.violinplot(data=dplot,**params['params_dist'],
    ax=ax,)
    _=[ax.text(**prms,ha='center') for prms in params['params_texts']]
    return ax


## FigureS10:panel0
def plot_hist_delta_interaction_score__between_species___within_Scer_(plotp="plot/hist_delta_interaction_score__between_species___within_Scer_.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from rohan.dandage.plot.annot import pval2annot
    from scipy.stats import norm
    ax=sns.distplot(dplot[params['colx']],bins=30,
    color=params['color'],
    fit_kws={"color":params['color'],'lw':2,'label':'fitted normal\ndistribution'},
    fit=norm,kde=False,ax=ax,)
    ax.axvline(dplot[params['colx']].mean(),
    color='red',linestyle='-',
    label=f"mean={dplot[params['colx']].mean():1.3f}")
    ax.axvline(dplot[params['colx']].median(),
    color='salmon',linestyle='-',
    label=f"median={dplot[params['colx']].median():1.3f}")    
    ax.set_xlabel(params['xlabel'])
    ax.set_ylabel('density')
    ax.set_xlim(-5,5)
    ax.legend(loc=1)
    return ax


## FigureS10:panel1
def plot_hist_delta_interaction_score__between_species___within_Suva_(plotp="plot/hist_delta_interaction_score__between_species___within_Suva_.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from rohan.dandage.plot.annot import pval2annot
    from scipy.stats import norm
    ax=sns.distplot(dplot[params['colx']],bins=30,
    color=params['color'],
    fit_kws={"color":params['color'],'lw':2,'label':'fitted normal\ndistribution'},
    fit=norm,kde=False,ax=ax,)
    ax.axvline(dplot[params['colx']].mean(),
    color='red',linestyle='-',
    label=f"mean={dplot[params['colx']].mean():1.3f}")
    ax.axvline(dplot[params['colx']].median(),
    color='salmon',linestyle='-',
    label=f"median={dplot[params['colx']].median():1.3f}")    
    ax.set_xlabel(params['xlabel'])
    ax.set_ylabel('density')
    ax.set_xlim(-5,5)
    ax.legend(loc=1)
    return ax


## Figure04:panel0
def plot_line_connections_similarity_interaction_types_pca(plotp="plot/line_connections_similarity_interaction_types_pca.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from matplotlib.lines import Line2D
    from rohan.dandage.plot.line import plot_connections
    ax=plot_connections(dplot=dplot,
    **{k:params[k] for k in params if k!='label2loc'},
    params_text={'ha':'center','va':'top'},
    ax=ax,test=False,)
    label2loc=params['label2loc']
    _=[ax.add_patch(patch) for patch in [plt.Rectangle(**label2loc[k]) for k in label2loc]]
    _=[ax.text(label2loc[k]['xy'][0],
    label2loc[k]['xy'][1]+label2loc[k]['height']*0.5,k,color=label2loc[k]['edgecolor'],va='center',ha='center') for k in ['parents','hybrid']]
    _=[ax.text(label2loc[k]['xy'][0]+label2loc[k]['width']*0.5,
    label2loc[k]['xy'][1]+label2loc[k]['height'],k,color=label2loc[k]['edgecolor'],va='bottom',ha='center') for k in ['intralogous PPIs']]
    return ax


## Figure04:panel1
def plot_network_hybrid_pca(plotp="plot/network_hybrid_pca.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot,ax=ax)
    plot_network_legend(ax,legend_size=False,**params['plot_network_legend'])
    return ax


## Figure04:panel2
def plot_scatter_interaction_score_Scer_intralogous_PPIs_in_hybrid_interaction_score_Suva_intralogous_PPIs_in_hybrid_interaction_score_interlogous_PPIs_pca(plotp="plot/scatter_interaction_score_Scer_intralogous_PPIs_in_hybrid_interaction_score_Suva_intralogous_PPIs_in_hybrid_interaction_score_interlogous_PPIs_pca.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()  
    dplot.plot.scatter(x=params['colx'],y=params['coly'],
    c=params['colz'],cmap='binary_r',
    ax=ax)
    _=ax.text(*ax.get_ylim(),params['label'],va='bottom')
    return ax


## Figure04:panel3
def plot_scatter_interaction_score_Scer_intralogous_PPIs_in_parent_interaction_score_Suva_intralogous_PPIs_in_parent_interaction_score_interlogous_PPIs_pca(plotp="plot/scatter_interaction_score_Scer_intralogous_PPIs_in_parent_interaction_score_Suva_intralogous_PPIs_in_parent_interaction_score_interlogous_PPIs_pca.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()  
    dplot.plot.scatter(x=params['colx'],y=params['coly'],
    c=params['colz'],cmap='binary_r',
    ax=ax)
    _=ax.text(*ax.get_ylim(),params['label'],va='bottom')
    return ax


## Figure04:panel4
def plot_scatter_interaction_score_interlogous_PPIs_predicted_pca(plotp="plot/scatter_interaction_score_interlogous_PPIs_predicted_pca.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from rohan.dandage.plot.ax_ import set_equallim
    from rohan.dandage.stat.corr import get_corr
    _=dplot.groupby('species type',sort=False).apply(lambda df: df.plot.scatter(x=params['colz'],
    y=f"{params['colz']} predicted",
    **{'color':params['element2color'][df.name],'label':'intralogous PPIs\nin '+df.name+'\n'+get_corr(df[params['colz']],
    df[f"{params['colz']} predicted"],
    method='spearman',
    bootstrapped=False,ci_type='max',outstr=True).replace('\n',' ')
    },
    ax=ax))
    ax.legend(title='predicted from',loc='center left',bbox_to_anchor=[1,0.5])
    ax.axis('scaled');ax.set(**{'xlim':[0,1],'ylim':[0,1]})
    ax=set_equallim(ax,diagonal=True)
    return ax


## FigureS12:panel0
def plot_network_parent_pca(plotp="plot/network_parent_pca.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot,ax=ax)
    plot_network_legend(ax,legend_size=False,**params['plot_network_legend'])
    return ax


## FigureS13:panel0
def plot_dist_delta_interaction_score_subset_divergence____(plotp="plot/dist_delta_interaction_score_subset_divergence____.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    from rohan.dandage.plot.dist import plot_dist_comparison
    if ax is None:ax=plt.subplot()
    ax.set_ylim(0,20)
    ax=plot_dist_comparison(df=dplot,**params,ax=ax)
    ax.grid(True)
    ax.set_xlim(-0.5,0.5)
    return ax


## Figure05:panel0
def plot_dist_predicted_from_intralogous_PPIs_in__predicted_actual__interaction_score_of_interlogous_PPIs_divergence_predicted_ms(plotp="plot/dist_predicted_from_intralogous_PPIs_in__predicted_actual__interaction_score_of_interlogous_PPIs_divergence_predicted_ms.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    from rohan.dandage.plot.dist import plot_dist_comparison
    if ax is None:ax=plt.subplot()
    ax=plot_dist_comparison(df=dplot,**params,ax=ax,params_legend={'bbox_to_anchor':[1,0.5],'loc':'center left'},
    kws_annot_boxplot={'yoff': -0.45},
    )
    ax.set_ylim(0,0.3)
    return ax


## Figure05:panel1
def plot_dist_protein_complex_description_interaction_score_ratio_zscore(plotp="plot/dist_protein_complex_description_interaction_score_ratio_zscore.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    from rohan.dandage.plot.ax_ import color_ticklabels
    if ax is None:ax=plt.subplot()
    sns.swarmplot(data=dplot,x=params['colx'],
    y='gene set description',hue='comparison type',palette=params['palette'],
    hue_order=params['hue_order'],
    ax=ax)
    ax.tick_params(axis='y', labelcolor=[0.25,0.25,0.25])
    ax=color_ticklabels(ax,params['yticklabel2color'],axis='y')
    ax.legend(loc='lower center',bbox_to_anchor=[0.5,1.01])
    ax.axvspan(2, 4, color=params['enrichmenttype2color']['high'], alpha=0.2,label='significantly high')
    ax.axvspan(-4, -2, color=params['enrichmenttype2color']['low'], alpha=0.2,label='significantly low')
    ax.set_xlim(-4,4)
    ax.set_ylabel('')
    ax.text(ax.get_xlim()[0],ax.get_ylim()[1],'low',color='b',ha='left',va='bottom')
    ax.text(ax.get_xlim()[1],ax.get_ylim()[1],'high',color='r',ha='right',va='bottom')    
    ax.set_xlabel('interaction score ratio\n(hybrid/parent)')
    return ax


## Figure05:panel2
def plot_dist_volcano_interaction_score_ratio_zscore_Biological_process(plotp="plot/dist_volcano_interaction_score_ratio_zscore_Biological_process.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from rohan.dandage.plot.scatter import plot_volcano
    ax=plot_volcano(dplot,**params,ax=ax,annots_off=0.4,annot_count_max=6)
    ax.text(ax.get_xlim()[0],ax.get_ylim()[1],'low',color='b',ha='left',va='bottom')
    ax.text(ax.get_xlim()[1],ax.get_ylim()[1],'high',color='r',ha='right',va='bottom')
    ax.set_xlabel('interaction score ratio\n(hybrid/parent)')    
    return ax


## Figure05:panel3
def plot_scatter_proteasome_parent_hybrid(plotp="plot/scatter_proteasome_parent_hybrid.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from rohan.dandage.plot.scatter import plot_scatter
    from rohan.dandage.plot.ax_ import set_equallim
    ax=dplot.plot.scatter(x=params['x'],y=params['y'],c=params['c'],
    cmap='RdBu_r',vmin=-2,vmax=2,
    ax=ax)
    ax=set_equallim(ax,diagonal=True)
    ax.grid(True)
    ax.text(params['annot']['x'],params['annot']['y'],params['annot']['s'],
    ha='center',va='bottom')
    return ax


## FigureS14:panel0
def plot_dist_paired_interlogous_PPIs_incompatibility(plotp="plot/dist_paired_interlogous_PPIs_incompatibility.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    dplot[params['coly']]=pd.Categorical(dplot[params['coly']],params['order'])
    from rohan.dandage.plot.dist import plot_dists
    ax=plot_dists(dplot,**params,
    xlims=None,
    cmap='Reds',
    palette=None,
    annot_pval=True,
    annot_n=True,
    annot_stat=False,
    params_dist={},
    params_violin={'scale':'count'},
    ax=ax)
    return ax


## FigureS15:panel0
def plot_dist_volcano_interaction_score_ratio_zscore_Molecular_function(plotp="plot/dist_volcano_interaction_score_ratio_zscore_Molecular_function.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from rohan.dandage.plot.scatter import plot_volcano
    ax=plot_volcano(dplot,**params,ax=ax,annots_off=0.4,annot_count_max=6)
    ax.text(ax.get_xlim()[0],ax.get_ylim()[1],'low',color='b',ha='left',va='bottom')
    ax.text(ax.get_xlim()[1],ax.get_ylim()[1],'high',color='r',ha='right',va='bottom')
    ax.set_xlabel('interaction score ratio\n(hybrid/parent)')    
    return ax


## FigureS15:panel1
def plot_dist_volcano_interaction_score_ratio_zscore_protein_complex(plotp="plot/dist_volcano_interaction_score_ratio_zscore_protein_complex.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    from rohan.dandage.plot.scatter import plot_volcano
    ax=plot_volcano(dplot,**params,ax=ax,annots_off=0.4,annot_count_max=6)
    ax.text(ax.get_xlim()[0],ax.get_ylim()[1],'low',color='b',ha='left',va='bottom')
    ax.text(ax.get_xlim()[1],ax.get_ylim()[1],'high',color='r',ha='right',va='bottom')
    ax.set_xlabel('interaction score ratio\n(hybrid/parent)')    
    return ax


## FigureS16:panel0
def plot_network_interactions_hybrid_CPX_554(plotp="plot/network_interactions_hybrid_CPX_554.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS16:panel1
def plot_network_interactions_parent_CPX_554(plotp="plot/network_interactions_parent_CPX_554.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS17:panel0
def plot_network_interactions_hybrid_CPX_1293(plotp="plot/network_interactions_hybrid_CPX_1293.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS17:panel1
def plot_network_interactions_parent_CPX_1293(plotp="plot/network_interactions_parent_CPX_1293.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS18:panel0
def plot_network_interactions_hybrid_CPX_3207(plotp="plot/network_interactions_hybrid_CPX_3207.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS18:panel1
def plot_network_interactions_parent_CPX_3207(plotp="plot/network_interactions_parent_CPX_3207.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS19:panel0
def plot_network_interactions_hybrid_CPX_1276(plotp="plot/network_interactions_hybrid_CPX_1276.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS19:panel1
def plot_network_interactions_parent_CPX_1276(plotp="plot/network_interactions_parent_CPX_1276.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS20:panel0
def plot_network_interactions_hybrid_CPX_1323(plotp="plot/network_interactions_hybrid_CPX_1323.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax


## FigureS20:panel1
def plot_network_interactions_parent_CPX_1323(plotp="plot/network_interactions_parent_CPX_1323.png",dplot=None,params=None,ax=None,fig=None,outd=None):
    plotp,dplot,params=get_plot_inputs(plotp=plotp,dplot=dplot,params=params,outd=f"{dirname(__file__)}");
    if ax is None:ax=plt.subplot()
    ax=plot_interactome_hybrid(dplot.reset_index(),title=None,test=False,ax=ax)
    ax.set_axis_off()
    return ax
