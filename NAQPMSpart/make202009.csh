#!/bin/csh
foreach MM (09)
  foreach DD (20 21 22 23 24 25 26 27 28 29 30 31)
  rm -rf  out/tmp/2020${MM}${DD}
  rm -rf  sm/tmp/2020${MM}${DD}
  rm -rf  drydep/tmp/2020${MM}${DD}
  rm -rf  wetdep/tmp/2020${MM}${DD}
  rm -rf  term/tmp/2020${MM}${DD}
  rm -rf  dust/tmp/2020${MM}${DD}
  rm -rf  seasalt/tmp/2020${MM}${DD}
  mkdir -p out/tmp/2020${MM}${DD}
  mkdir -p sm/tmp/2020${MM}${DD}
  mkdir -p drydep/tmp/2020${MM}${DD}
  mkdir -p wetdep/tmp/2020${MM}${DD}
  mkdir -p term/tmp/2020${MM}${DD}
  mkdir -p dust/tmp/2020${MM}${DD}
  mkdir -p seasalt/tmp/2020${MM}${DD}
  end
end
foreach MM (10)
  foreach DD (01 02 03 04 05 06 07 08 09 10)
  rm -rf  out/tmp/2020${MM}${DD}
  rm -rf  sm/tmp/2020${MM}${DD}
  rm -rf  drydep/tmp/2020${MM}${DD}
  rm -rf  wetdep/tmp/2020${MM}${DD}
  rm -rf  term/tmp/2020${MM}${DD}
  rm -rf  dust/tmp/2020${MM}${DD}
  rm -rf  seasalt/tmp/2020${MM}${DD}
  mkdir -p out/tmp/2020${MM}${DD}
  mkdir -p sm/tmp/2020${MM}${DD}
  mkdir -p drydep/tmp/2020${MM}${DD}
  mkdir -p wetdep/tmp/2020${MM}${DD}
  mkdir -p term/tmp/2020${MM}${DD}
  mkdir -p dust/tmp/2020${MM}${DD}
  mkdir -p seasalt/tmp/2020${MM}${DD}
  end
end

