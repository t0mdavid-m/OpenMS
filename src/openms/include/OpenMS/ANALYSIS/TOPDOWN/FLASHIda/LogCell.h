// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <string>

namespace OpenMS
{
  /**
   * @brief The cell grammar of FLASHIda's four TSV logs: one writer rule and its inverse.
   *
   * IdaLogger writes every cell that carries a tracking id -- the ids, the space-joined id lists and
   * scan_description -- through quoteLogCell(); no other cell can begin with '"' (numbers, fixed
   * vocabularies, ion labels, ProForma). The readers of those logs in this repository,
   * ScanCommandJoin::parse (FLASHDeconv's -FD:scan_commands) and the class tests' TSVFile, read them
   * back through unquoteLogId() / unquoteLogCell().
   *
   * Why a grammar at all: tracking ids are base-94 over 0x21-0x7E (ScanCommandQueue::encode), and
   * that alphabet includes '"'. Every id from 8,836 to 17,671 BEGINS with it -- over a third of the
   * commands of a 14,000-command run -- and a CSV-style reader (Python's csv and pandas at their
   * defaults) takes a '"' at the start of a cell as an opening quote: it swallows tabs and line ends up
   * to the next '"' and shifts or merges rows, silently or with a parse error.
   *
   * The rule is RFC 4180's, applied where those readers act on a quote: a cell that BEGINS with '"'
   * is enclosed in quotes with each of its quotes doubled, so "!# is written """!#". Every other
   * cell is written as it is: a '"' anywhere else is read literally by those readers, and leaving
   * those cells alone keeps every log of a run below 8,836 commands byte-identical. No cell carries a
   * tab or a line end (the id alphabet and the scan_description grammar exclude both), so nothing else
   * needs quoting and splitting a row on tabs stays exact; a tab-splitting reader only has to unquote.
   */
  inline std::string quoteLogCell(const std::string& cell)
  {
    if (cell.empty() || cell.front() != '"') { return cell; }
    std::string out;
    out.reserve(cell.size() + 4);
    out.push_back('"');
    for (const char c : cell)
    {
      out.push_back(c);
      if (c == '"') { out.push_back('"'); }
    }
    out.push_back('"');
    return out;
  }

  /**
   * @brief Inverse of quoteLogCell(), and of RFC 4180 quoting in general.
   *
   * A cell enclosed in quotes whose inner quotes are all doubled is unquoted; anything else is
   * returned unchanged. The general form matters because a log re-saved by pandas or Excel quotes
   * every cell that CONTAINS a '"' (!!" becomes "!!"""), not only the ones quoteLogCell() quotes.
   *
   * Exact on any file written with quoteLogCell(), where no unquoted cell begins with '"'. A file
   * written BEFORE the rule holds raw 3-character ids such as "!" (id 8,837), which this function
   * reads as a quoted !; read tracking ids with unquoteLogId(), which is exact on both.
   */
  inline std::string unquoteLogCell(const std::string& cell)
  {
    if (cell.size() < 2 || cell.front() != '"' || cell.back() != '"') { return cell; }
    std::string out;
    out.reserve(cell.size() - 2);
    for (std::size_t i = 1; i + 1 < cell.size(); ++i)   // the interior, [1, size - 2]
    {
      if (cell[i] == '"')
      {
        if (i + 2 < cell.size() && cell[i + 1] == '"') { out.push_back('"'); ++i; continue; }
        return cell;   // a lone inner quote: not a quoted cell
      }
      out.push_back(cell[i]);
    }
    return out;
  }

  /// unquoteLogCell() for a tracking-id cell, exact for files written before AND after the rule: a raw
  /// id is always 3 characters (ScanCommandQueue::encode), a quoted one at least 6 ("!! -> """!!"),
  /// so a 3-character cell is never unquoted -- not even "!", which is the raw id 8,837.
  inline std::string unquoteLogId(const std::string& cell)
  {
    return cell.size() <= 3 ? cell : unquoteLogCell(cell);
  }
} // namespace OpenMS
