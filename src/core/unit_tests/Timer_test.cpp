/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <unistd.h>

#define BOOST_TEST_MODULE Utils::Timing::Timer test
#include <boost/test/included/unit_test.hpp>

#include "../utils/Timer.hpp"

using namespace Utils::Timing;
using std::cout;
using std::endl;

BOOST_AUTO_TEST_CASE(sampling) {
  Timer &t = Timer::get_timer("test_timer");

  for (int i = 0; i < 100; i++) {
    t.start();
    usleep(1000);
    t.stop();
  }

  BOOST_CHECK(t.stats().n() == 100);
}

/**
 * Check that we can get the timer back by name
 */
BOOST_AUTO_TEST_CASE(get_timer) {
  Timer &t = Timer::get_timer("test_timer");

  t.start();
  usleep(1000);
  t.stop();

  BOOST_CHECK(t.stats().n() == 101);
}

BOOST_AUTO_TEST_CASE(reset_timer) {
  Timer &t = Timer::get_timer("test_timer");
  BOOST_CHECK(t.stats().n() == 101);
  t.reset();
  BOOST_CHECK(t.stats().n() == 0);
}

BOOST_AUTO_TEST_CASE(reset_all_timer) {
  Timer &t = Timer::get_timer("test_timer");
  Timer &u = Timer::get_timer("test_timer_2");

  t.start();
  t.stop();

  u.start();
  u.stop();

  t.start();
  t.stop();

  t.start();
  t.stop();

  Timer::reset_all();

  BOOST_CHECK((Timer::get_timer("test_timer").stats().n() == 0) &&
              (Timer::get_timer("test_timer_2").stats().n() == 0));
}
