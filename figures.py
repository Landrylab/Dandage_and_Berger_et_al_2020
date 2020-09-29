from rohan.global_imports import *
from rohan.dandage.figs.figure import *
from rohan.dandage.plot.schem import *
from .plots import *
import warnings
warnings.filterwarnings("ignore")

def Figure01(ind,outd):
    fig=plt.figure(figsize=[8,10],linewidth=5)
    axs=[
         plt.subplot2grid([1,5],[0,0],1,3),
         plt.subplot2grid([1,5],[0,3],1,2),
    ]
    ps=['plot/schem_protein_complex_possibilities_of_assembly_scenarios.svg',
        'plot/schem01_method.svg']
    for p,ax in zip(ps,axs):
        plot_schem(p,ax=ax,force=False)
        ax.axis('off')
    plt.subplots_adjust(wspace=0.3, 
                       )
    labelplots(fig=fig,axes=axs,params_alignment={(0,1):['y']})
    savefig(f'{outd}/figures/Figure01',tight_layout=True,fmts=['png','svg'])

def FigureS01(ind,outd):
    fig,axs=plt.subplots(nrows=2,ncols=2,figsize=[10,5.5],linewidth=5)
    axs=np.ravel(axs)
    plot_hist_orthologs_protein_sequence_divergence____(ax=axs[0])
    plot_hist_orthologs_dissimilar_peptides_per_protein__Jaccard_distance___(ax=axs[2])
    plot_bar_peptides_detected_by_unique_not_unique(ax=axs[1])
    plot_bar_proteins_detected_by_unique_not_unique(ax=axs[3])
    labelsubplots(axes=[axs[0],axs[2],axs[1],axs[3]],
                  xoff=-0.01,yoff=0.05)
    plt.tight_layout()
    savefig(f'{outd}/figures/FigureS01',tight_layout=True,fmts=['png','svg'])

def FigureS02(ind,outd):
    plt.figure(figsize=[8,8],linewidth=5)
    ax=plt.subplot()
    plot_schem(f'{ind}/plot/schem02_method_massspec.svg',ax=ax,force=False)
    savefig(f'{outd}/figures/FigureS02',tight_layout=True,fmts=['png','svg'])

def FigureS03(ind,outd):
    plt.figure(figsize=[5,5])
    ax=plot_heatmap_interaction_qc_interaction_bool_db_ROC_AUC()
    savefig(f'{outd}/figures/FigureS03',tight_layout=True,fmts=['png','svg'])

def Figure02(ind,outd):
    def plot_line_interaction_score_protein_dosage_balance_(ax):
        ## panel
        funs=[
        plot_line_interaction_score_protein_dosage_balance_within_parent_Scer_Scer,
        plot_line_interaction_score_protein_dosage_balance_within_parent_Suva_Suva,
        plot_line_interaction_score_protein_dosage_balance_within_hybrid_Scer_Scer,
        plot_line_interaction_score_protein_dosage_balance_within_hybrid_Suva_Suva,
        plot_line_interaction_score_protein_dosage_balance_between_hybrid_Scer_Suva,
          ]
        [f(ax=ax) for f in funs]
        ax.legend(bbox_to_anchor=[1,1])
        ax.set_xlim(0.5,0.96)
        return ax   
    plt.figure(figsize=[9,10.5],linewidth=5)
    axs=[
         plt.subplot2grid([3,2],[0,0],1,1),
         plt.subplot2grid([3,2],[0,1],1,1),
         plt.subplot2grid([3,2],[1,0],1,1),
         plt.subplot2grid([3,2],[1,1],1,1),
         plt.subplot2grid([3,2],[2,0],1,1),
         plt.subplot2grid([3,2],[2,1],1,1),
        ]
    funs=[
        plot_line_interaction_qc_ROC_AUC_DTW__window_04__score__quantile_0_75__interaction_bool_db_formated,
        plot_scatter_protein_abundance_db,
        plot_violin_interaction_score_mean_by_protein_complexes_and_random,
        plot_violin_interaction_score_mean_by_protein_subcomplex,
        plot_scatter_interaction_score_distance_proteasome,
        plot_line_interaction_score_protein_dosage_balance_,
    ]
    for fun,ax in zip(funs,axs):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            ax=fun(ax=ax)
    labelsubplots(axes=axs,xoff=-0.01)
    plt.tight_layout()
    savefig(f'{outd}/figures/Figure02',tight_layout=True,fmts=['png','svg'])

def FigureS04(ind,outd):
    plt.figure(figsize=[7,6])
    plot_network_within_parent_Scer_Scer_proteasome_CPX_2262()
    savefig(f'{outd}/figures/FigureS04',tight_layout=True,fmts=['png','svg'])

def FigureS05(ind,outd):
    plt.figure(figsize=[3,7])
    axs=[
        plt.subplot2grid([2,1],[0,0],1,1),
        plt.subplot2grid([2,1],[1,0],1,1),
        ]
    funs=[
    plot_line_interaction_score_protein_complexes_and_random,
    plot_line_protein_dosage_balance_protein_complexes_and_random,
    ]
    for fun,ax in zip(funs,np.ravel(axs)):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            ax=fun(ax=ax)
    labelsubplots(axes=np.ravel(axs),xoff=0)
    plt.subplots_adjust(
        wspace=0.1, hspace=0.3)
    savefig(f'{outd}/figures/FigureS05',tight_layout=True,fmts=['png','svg'])

def Figure03(ind,outd):
    fig=plt.figure(figsize=[13,18],linewidth=5)
    axs=[
         plt.subplot2grid([9,5],[0,0],3,3,),
         plt.subplot2grid([9,5],[0,3],3,2,),
         plt.subplot2grid([9,3],[3,0],2,1,),
         plt.subplot2grid([9,3],[3,1],2,1,),
         plt.subplot2grid([9,3],[3,2],2,1,),
         plt.subplot2grid([9,4],[5,1],3,2,),
        ]
    def plot_schem01(ax): return plot_schem(f'{ind}/plot/schem_additivity_between_species_interactions.svg',ax=ax,force=False)
    funs=[
        plot_line_connections_similarity_interaction_types_ms,
        plot_schem01,
        plot_contour_within_hybrids__Scer__Scer__within_hybrids__Suva__Suva__between_hybrids__Scer__Suva_,
        plot_contour_within_parents__Scer__Scer__within_parents__Suva__Suva__between_hybrids__Scer__Suva_,  
        plot_scatter_interaction_score_interlogous_PPIs_predicted_ms,
        plot_circle_complexes_hybrid,
    ]
    for fun,ax in zip(funs,axs):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            if not 'contour' in fun.__name__:
                ax=fun(ax=ax)
            else:
                ax=fun(ax=ax,fig=fig,params={'params_cbar':{'bbox_to_anchor':(1.02, -0.8, 0.8, 1.75)}})            
    labelsubplots(axes=np.ravel(axs),xoff=0,yoff=0.02)
    plt.subplots_adjust(
        wspace=0.75, hspace=2.5)
    savefig(f'{outd}/figures/Figure03',tight_layout=True,fmts=['png','svg'])

def FigureS06(ind,outd):
    from rohan.dandage.plot.colors import reset_legend_colors
    funs=[
    plot_dendogram_interaction_score_within_species,
    plot_dendogram_protein_abundance,
    plot_heatmap_interaction_score_within_species,
    plot_heatmap_protein_abundance,
    ]
    fig,axes=plt.subplots(ncols=2,nrows=2,figsize=[6,6],
                          gridspec_kw={'height_ratios': [1, 3]},
                         linewidth=5)
    for fun,ax in zip(funs,np.ravel(axes)):
        ax=fun(ax=ax)
        if 'heatmap' in fun.__name__:
            ax.set_xlabel('')
    labelsubplots(axes=np.ravel(axes)[:2],yoff=-0.1,xoff=0.1)
    plt.subplots_adjust(
        left=0.075, right=1.075, 
        bottom=0.075,top=0.95,  
        wspace=0.5, hspace=0.5)
    savefig(f'{outd}/figures/FigureS06',tight_layout=True,fmts=['png','svg'])

def FigureS07(ind,outd):
    plt.figure(figsize=[10,15],linewidth=5)
    axs=[
        plt.subplot2grid([2,1],[0,0],1,2),
        plt.subplot2grid([2,1],[1,0],1,2),
        ]
    funs=[plot_network_interactions_hybrid_CPX_2262,
    plot_network_interactions_parent_CPX_2262,
         ]
    for fun,ax in zip(funs,axs):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            ax=fun(ax=ax)
            if 'network' in fun.__name__:
                ax.set_title(fun.__name__.split('_')[-3],loc='left')
    labelsubplots(axes=np.ravel(axs),xoff=0)
    plot_network_legend(axs[1])
    savefig(f'{outd}/figures/FigureS07',tight_layout=True,fmts=['png','svg'])

def FigureS08(ind,outd):
    plt.figure(figsize=[8,11])
    axs=[plt.subplot2grid([2,1],[0,0],1,1),
         plt.subplot2grid([2,1],[1,0],1,1)
        ]
    funs=[plot_network_interactions_hybrid_CPX_1671,
    plot_network_interactions_parent_CPX_1671]
    for fun,ax in zip(funs,axs):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            ax=fun(ax=ax)
            ax.set_title(fun.__name__.split('_')[-3],loc='left')
    labelsubplots(axes=np.ravel(axs),xoff=0)
    plot_network_legend(axs[1])
    plt.tight_layout()
    savefig(f'{outd}/figures/FigureS08',tight_layout=True,fmts=['png','svg'])

def FigureS09(ind,outd):
    plt.figure(figsize=[3,3])
    plot_dist_interaction_type_interaction_score_complex_counts()
    savefig(f'{outd}/figures/FigureS09',tight_layout=True,fmts=['png','svg'])

def FigureS10(ind,outd):
    plt.figure(figsize=[5,6],linewidth=5)
    axs=[
         plt.subplot2grid([2,1],[0,0],1,1,),
         plt.subplot2grid([2,1],[1,0],1,1,),]
    funs=[plot_hist_delta_interaction_score__between_species___within_Scer_,
          plot_hist_delta_interaction_score__between_species___within_Suva_,
         ]
    for fun,ax in zip(funs,axs):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            ax=fun(ax=ax)
    labelsubplots(axes=np.ravel(axs),xoff=0)
    plt.tight_layout()
    savefig(f'{outd}/figures/FigureS10',tight_layout=True,fmts=['png','svg'])

def Figure04(ind,outd):
    fig=plt.figure(figsize=[13,16])
    axs=[
        plt.subplot2grid([16,5],[0,0],5,2,),
        plt.subplot2grid([16,5],[0,2],4,3,),
        plt.subplot2grid([16,6],[6,0],5,3,),
        plt.subplot2grid([16,12],[7,7],3,3,),
        plt.subplot2grid([16,12],[12,1],3,3,),
        plt.subplot2grid([16,6],[12,3],3,2,),
    ]
    def plot_schem01(ax):return plot_schem(f'{ind}/plot/schem_pca_summary.svg',margin=0,ax=ax)
    funs=[
    plot_schem01,
    plot_network_hybrid_pca,
    plot_line_connections_similarity_interaction_types_pca,
    plot_scatter_interaction_score_Scer_intralogous_PPIs_in_hybrid_interaction_score_Suva_intralogous_PPIs_in_hybrid_interaction_score_interlogous_PPIs_pca,
    plot_scatter_interaction_score_Scer_intralogous_PPIs_in_parent_interaction_score_Suva_intralogous_PPIs_in_parent_interaction_score_interlogous_PPIs_pca,
    plot_scatter_interaction_score_interlogous_PPIs_predicted_pca,    
    ]
    for fun,ax in zip(funs,np.ravel(axs)):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            ax=fun(ax=ax)
    labelplots(fig=fig,axes=np.ravel(axs),params_alignment={(0,1):['y']},yoff=0.01)
    plt.subplots_adjust(
        hspace=0.5)
    savefig(f'{outd}/figures/Figure04',tight_layout=True,fmts=['png','svg'])

def FigureS11(ind,outd):
    plt.figure(figsize=[8,8])
    ax=plt.subplot()
    plot_schem(f'{ind}/plot/schem_pca_detailed.svg',ax=ax)
    savefig(f'{outd}/figures/FigureS11',tight_layout=True,fmts=['png','svg'])

def FigureS12(ind,outd):
    plot_network_parent_pca()
    savefig(f'{outd}/figures/FigureS12',tight_layout=True,fmts=['png','svg'])

def FigureS13(ind,outd):
    plt.figure(figsize=[11,7],linewidth=5)
    axs=[
         plt.subplot2grid([3,1],[0,0],2,1,),
         plt.subplot2grid([3,3],[2,1],1,1,),
        ]
    def plot_schema(ax):
        return plot_schem(f'{ind}/plot/schem_delta_interaction_ortholog_seq_identity.svg',ax=ax,force=False)
    funs=[plot_schema,
         plot_dist_delta_interaction_score_subset_divergence____]
    for fun,ax in zip(funs,axs):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            ax=fun(ax=ax)
    labelsubplots(axes=np.ravel(axs),xoff=0.1)
    plt.subplots_adjust(
        hspace=0.4)
    savefig(f'{outd}/figures/FigureS13',tight_layout=True,fmts=['png','svg'])

def Figure05(ind,outd):
    plt.figure(figsize=[16,15])
    axs=[
         plt.subplot2grid([7,20],[0,3],2,7,),
         plt.subplot2grid([14,20],[4,5],5,4,),
         plt.subplot2grid([20,5],[1,4],19,1),
         plt.subplot2grid([14,20],[10,4],4,6),
        ]
    funs=[
        plot_dist_predicted_from_intralogous_PPIs_in__predicted_actual__interaction_score_of_interlogous_PPIs_divergence_predicted_ms,
        plot_dist_volcano_interaction_score_ratio_zscore_Biological_process,   
        plot_dist_protein_complex_description_interaction_score_ratio_zscore,
        plot_scatter_proteasome_parent_hybrid,
         ]
    for fun,ax in zip(funs,axs):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            ax=fun(ax=ax)
    labelsubplots(axes=np.ravel(axs),xoff=0)
    plt.subplots_adjust(wspace=0.2, hspace=1.5)
    savefig(f'{outd}/figures/Figure05',tight_layout=True,fmts=['png','svg'])

def FigureS14(ind,outd):
    plt.figure(figsize=[3,3])
    ax=plt.subplot()
    plot_dist_paired_interlogous_PPIs_incompatibility(ax=ax)
    savefig(f'{outd}/figures/FigureS14',tight_layout=True,fmts=['png','svg'])

def FigureS15(ind,outd):
    fig,axs=plt.subplots(nrows=2,ncols=1,figsize=[3,12],sharex=False)
    funs=[
          plot_dist_volcano_interaction_score_ratio_zscore_protein_complex,
        plot_dist_volcano_interaction_score_ratio_zscore_Molecular_function,
         ]
    for fun,ax in zip(funs,axs):
        if isinstance(fun,list):
            [f(ax=ax) for f in fun]
        else:
            ax=fun(ax=ax)
    labelsubplots(axes=np.ravel(axs),xoff=0)
    plt.subplots_adjust(
        wspace=0, hspace=0.3)
    savefig(f'{outd}/figures/FigureS15',tight_layout=True,fmts=['png','svg'])

def FigureS16(ind,outd):
    plot_complex_hybrid_parent([globals()[k] for k in ['plot_network_interactions_parent_CPX_554',
                               'plot_network_interactions_hybrid_CPX_554']],
                              figsize=[7,9])
    savefig(f'{outd}/figures/FigureS16',tight_layout=True,fmts=['png','svg'])

def FigureS17(ind,outd):
    plot_complex_hybrid_parent([globals()[k] for k in ['plot_network_interactions_parent_CPX_1293',
                                'plot_network_interactions_hybrid_CPX_1293'
                               ]],figsize=[7,9])
    savefig(f'{outd}/figures/FigureS17',tight_layout=True,fmts=['png','svg'])

def FigureS18(ind,outd):
    plot_complex_hybrid_parent([globals()[k] for k in ['plot_network_interactions_parent_CPX_3207',
                                'plot_network_interactions_hybrid_CPX_3207',
                               ]],
                              figsize=[7,9])
    savefig(f'{outd}/figures/FigureS18',tight_layout=True,fmts=['png','svg'])

def FigureS19(ind,outd):
    plot_complex_hybrid_parent([globals()[k] for k in ['plot_network_interactions_parent_CPX_1276',
                               'plot_network_interactions_hybrid_CPX_1276']],
                              figsize=[7,9])
    savefig(f'{outd}/figures/FigureS19',tight_layout=True,fmts=['png','svg'])

def FigureS20(ind,outd):
    plot_complex_hybrid_parent([globals()[k] for k in ['plot_network_interactions_parent_CPX_1323',
                               'plot_network_interactions_hybrid_CPX_1323']],
                              figsize=[7,9])
    savefig(f'{outd}/figures/FigureS20',tight_layout=True,fmts=['png','svg'])