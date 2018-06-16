import re
RComp=re.compile
runID = RComp('[\^/](\d{1,4}[RC][_/].*)')
baseSpace = '(?P<lead>.*)_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)'