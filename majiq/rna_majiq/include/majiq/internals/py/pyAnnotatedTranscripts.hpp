/**
 * pyAnnotatedTranscripts.hpp
 *
 * Copyright 2025 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYANNOTATEDTRANSCRIPTS_HPP
#define MAJIQ_PYBIND_PYANNOTATEDTRANSCRIPTS_HPP

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <memory>
#include <utility>
#include <vector>
#include "../Exons.hpp"
#include "../GeneIntrons.hpp"
#include "coordinatesProperties.hpp"
#include "templateTypes.hpp"
#include "../AnnotatedTranscripts.hpp"

namespace majiq {
namespace bindings {

using pyAnnotatedTranscripts_t = pyClassShared_t<majiq::AnnotatedTranscripts>;

inline void init_AnnotatedTranscripts(pyAnnotatedTranscripts_t& pyAnnotatedTranscripts) {
    define_coordinates_properties(pyAnnotatedTranscripts);
    pyAnnotatedTranscripts
        .def(pybind11::init([](std::shared_ptr<majiq::Genes> genes,
                               pybind11::array_t<size_t> _gene_idx,
                               pybind11::array_t<position_t> _start,
                               pybind11::array_t<position_t> _end,
                               pybind11::list transcript_id,
                               pybind11::array_t<size_t> _exon_idx_start,
                               pybind11::array_t<size_t> _exon_idx_end,
                               pybind11::array_t<position_t> _exon_start,
                               pybind11::array_t<position_t> _exon_end
                               ) {
               if (_gene_idx.ndim() != 1 || _start.ndim() != 1 ||
                   _end.ndim() != 1) {
                 throw std::invalid_argument(
                     "AnnotatedTranscripts::init input arrays must be 1D");
               }
               if (_gene_idx.shape(0) != _start.shape(0) ||
                   _gene_idx.shape(0) != _end.shape(0)) {
                 throw std::invalid_argument(
                     "AnnotatedTranscripts::init input arrays have conflicting sizes");
               }
               // unchecked accesses to numpy array
               auto gene_idx = _gene_idx.unchecked<1>();
               auto start = _start.unchecked<1>();
               auto end = _end.unchecked<1>();

               auto exon_idx_start = _exon_idx_start.cast<std::vector<size_t>>();
               auto exon_idx_end = _exon_idx_end.cast<std::vector<size_t>>();
               auto exon_start = _exon_start.cast<std::vector<position_t>>();
               auto exon_end = _exon_end.cast<std::vector<position_t>>();

               // create vector of genes matching input arrays
               std::vector<majiq::AnnotatedTranscript> transcript_vec{};
               transcript_vec.reserve(gene_idx.shape(0));
               for (pybind11::ssize_t i = 0; i < gene_idx.shape(0); ++i) {
                 transcript_vec.push_back(
                     majiq::AnnotatedTranscript{
                         majiq::KnownGene{gene_idx(i), genes},
                         majiq::ExonIntervalT{start(i), end(i)},
                         transcript_id[i].cast<majiq::geneid_t>(),
                         std::vector<ClosedInterval>()
                     });
               }
               pybind11::gil_scoped_release release;  // release GIL at this stage
               return std::make_shared<AnnotatedTranscripts>(genes, std::move(transcript_vec),
                   exon_idx_start, exon_idx_end, exon_start, exon_end);
             }),
             "Create AnnotatedTranscripts object using Genes and info about each transcript",
             pybind11::arg("genes"), pybind11::arg("gene_idx"),
             pybind11::arg("start"), pybind11::arg("end"),
             pybind11::arg("transcript_id"),
             pybind11::arg("exon_idx_start"), pybind11::arg("exon_idx_end"),
             pybind11::arg("exon_start"), pybind11::arg("exon_end"))
    .def_property_readonly(
          "transcript_id",
          &AnnotatedTranscripts::transcript_ids,
          pybind11::call_guard<pybind11::gil_scoped_release>(),
          "array[str] of annotated transcript ids")
    .def_property_readonly(
        "exon_idx_start",
        [](pybind11::object& annotatedTranscripts_object) -> pybind11::array_t<size_t> {
          AnnotatedTranscripts& annotated_transcripts = annotatedTranscripts_object.cast<AnnotatedTranscripts&>();
          return pybind11::array(annotated_transcripts.exon_idx_start.size(), annotated_transcripts.exon_idx_start.data());
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[int] of indexes indicating beginning of exon_start/exon_end arrays for this transcript")
    .def_property_readonly(
        "exon_idx_end",
        [](pybind11::object& annotatedTranscripts_object) -> pybind11::array_t<size_t> {
          AnnotatedTranscripts& annotated_transcripts = annotatedTranscripts_object.cast<AnnotatedTranscripts&>();
          return pybind11::array(annotated_transcripts.exon_idx_end.size(), annotated_transcripts.exon_idx_end.data());
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[int] of indexes indicating end of exon_start/exon_end arrays for this transcript")
    .def_property_readonly(
        "exon_start",
        [](pybind11::object& annotatedTranscripts_object) -> pybind11::array_t<position_t> {
          AnnotatedTranscripts& annotated_transcripts = annotatedTranscripts_object.cast<AnnotatedTranscripts&>();
          return pybind11::array(annotated_transcripts.exon_start.size(), annotated_transcripts.exon_start.data());
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[int] of annotated exon starts")
    .def_property_readonly(
        "exon_end",
        [](pybind11::object& annotatedTranscripts_object) -> pybind11::array_t<position_t> {
          AnnotatedTranscripts& annotated_transcripts = annotatedTranscripts_object.cast<AnnotatedTranscripts&>();
          return pybind11::array(annotated_transcripts.exon_end.size(), annotated_transcripts.exon_end.data());
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[int] of annotated exon ends")
    ;
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYANNOTATEDTRANSCRIPTS_HPP
