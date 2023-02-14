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
#pragma once

#include <algorithm>
#include <cassert>
#include <cstdarg>
#include <fstream>
#include <iostream>

namespace music {

enum log_level : int {
  off     = 0,
  fatal   = 1,
  error   = 2,
  warning = 3,
  info    = 4,
  user    = 5,
  debug   = 6,
};

class logger {
private:
  static log_level log_level_;
  static std::ofstream output_file_;

public:
  logger()  = default;
  ~logger() = default;

  static void set_level(const log_level &level);
  static log_level get_level();

  static void set_output(const std::string filename);
  static void unset_output();

  static std::ofstream &get_output();

  template <typename T> logger &operator<<(const T &item) {
    std::cout << item;
    if (output_file_.is_open()) {
      output_file_ << item;
    }
    return *this;
  }

  logger &operator<<(std::ostream &(*fp)(std::ostream &)) {
    std::cout << fp;
    if (output_file_.is_open()) {
      output_file_ << fp;
    }
    return *this;
  }
};

class log_stream {
private:
  logger &logger_;
  log_level stream_level_;
  std::string line_prefix_, line_postfix_;

  bool newline;

public:
  log_stream(logger &logger, const log_level &level)
    : logger_(logger), stream_level_(level), newline(true) {
    switch (stream_level_) {
      case log_level::fatal:
        line_prefix_ = "\033[31mFatal : ";
        break;
      case log_level::error:
        line_prefix_ = "\033[31mError : ";
        break;
      case log_level::warning:
        line_prefix_ = "\033[33mWarning : ";
        break;
      case log_level::info:
      case log_level::user:
        //line_prefix_ = " | Info    | ";
        line_prefix_ = " \033[0m";
        break;
      case log_level::debug:
        line_prefix_ = "Debug : \033[0m";
        break;
      default:
        line_prefix_ = "\033[0m";
        break;
    }
    line_postfix_ = "\033[0m";
  }
  ~log_stream() = default;

  inline std::string GetPrefix() const {
    return line_prefix_;
  }

  template <typename T> log_stream &operator<<(const T &item) {
    if (logger::get_level() >= stream_level_) {
      if (newline) {
        logger_ << line_prefix_;
        newline = false;
      }
      logger_ << item;
    }
    return *this;
  }

  log_stream &operator<<(std::ostream &(*fp)(std::ostream &)) {
    if (logger::get_level() >= stream_level_) {
      logger_ << fp;
      logger_ << line_postfix_;
      newline = true;
    }
    return *this;
  }

  inline void Print(const char *str, ...) {
    char out[1024];
    va_list argptr;
    va_start(argptr, str);
    vsprintf(out, str, argptr);
    va_end(argptr);
    std::string out_string = std::string(out);
    out_string.erase(std::remove(out_string.begin(), out_string.end(), '\n'),
                     out_string.end());
    (*this) << out_string << std::endl;
  }
};

// global instantiations for different levels
extern logger glogger;
extern log_stream flog;
extern log_stream elog;
extern log_stream wlog;
extern log_stream ilog;
extern log_stream ulog;
extern log_stream dlog;

} // namespace music
