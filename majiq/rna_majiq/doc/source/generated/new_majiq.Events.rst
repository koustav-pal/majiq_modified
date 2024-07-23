new\_majiq.Events
=================

.. currentmodule:: rna_majiq

.. autoclass:: Events

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~Events.__init__
      ~Events.broadcast_eidx_to_ecidx
      ~Events.connection_contig_idx
      ~Events.connection_denovo
      ~Events.connection_end
      ~Events.connection_gene_idx
      ~Events.connection_other_exon_idx
      ~Events.connection_start
      ~Events.connections_slice_for_event
      ~Events.e_idx_slice_for_gene
      ~Events.ec_dataframe
      ~Events.ec_idx_slice_for_event
      ~Events.event_has_alt_exons
      ~Events.event_has_other_alt_ss
      ~Events.event_has_ref_alt_ss
      ~Events.event_id
      ~Events.event_legacy_a3ss
      ~Events.event_legacy_a5ss
      ~Events.from_arrays
      ~Events.from_zarr
      ~Events.has_intron
      ~Events.index
      ~Events.merge_dataframes
      ~Events.select_eidx_to_select_ecidx
      ~Events.slice_for_gene
      ~Events.to_zarr
      ~Events.unique_events_mask
   
   

   
   
   .. rubric:: Attributes

   .. autosummary::
   
      ~Events.connection_e_idx
      ~Events.connection_idx
      ~Events.connection_is_intron
      ~Events.connection_ref_exon_idx
      ~Events.contigs
      ~Events.df
      ~Events.df_event_connections
      ~Events.df_events
      ~Events.e_idx
      ~Events.e_idx_end
      ~Events.e_idx_start
      ~Events.ec_idx
      ~Events.ec_idx_end
      ~Events.ec_idx_start
      ~Events.event_size
      ~Events.event_type
      ~Events.exons
      ~Events.genes
      ~Events.introns
      ~Events.is_intron
      ~Events.junctions
      ~Events.num_connections
      ~Events.num_events
      ~Events.ref_exon_idx
      ~Events.save_df
   
   