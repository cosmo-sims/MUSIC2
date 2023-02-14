// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn & Michael Michaux (this file)
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <logger.hh>

namespace music {

std::ofstream logger::output_file_;
log_level logger::log_level_ = log_level::off;

void logger::set_level(const log_level &level) {
  log_level_ = level;
}

log_level logger::get_level() {
  return log_level_;
}

void logger::set_output(const std::string filename) {
  if (output_file_.is_open()) {
    output_file_.close();
  }
  output_file_.open(filename, std::ofstream::out);
  assert(output_file_.is_open());
}

void logger::unset_output() {
  if (output_file_.is_open()) {
    output_file_.close();
  }
}

std::ofstream &logger::get_output() {
  return output_file_;
}

// global instantiations for different levels
logger the_logger;
log_stream flog(the_logger, log_level::fatal);
log_stream elog(the_logger, log_level::error);
log_stream wlog(the_logger, log_level::warning);
log_stream ilog(the_logger, log_level::info);
log_stream ulog(the_logger, log_level::user);
log_stream dlog(the_logger, log_level::debug);

} // namespace music
