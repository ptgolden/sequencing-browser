// jshint laxcomma: true, asi: true

/* Constants */

var consts = {
  // Height of total plot area
  PLOT_WIDTH: 900,
  PLOT_HEIGHT: 610,

  // Height and width of the cell graphs
  CELL_HEIGHT: 500,
  CELL_WIDTH: 140,

  // Padding above and below cell graphs
  CELL_PADDING: 30
};

var DATA_FILES = [
  'data/193_4cell_RPKM.txt',
  'data/194_4cell_RPKM.txt',
  'data/195_4cell_RPKM.txt',
  'data/196_4cell_RPKM.txt'
];


/* Utilities */

function parseFile(dataFile) {
  return new Promise(function (resolve, reject) {
    d3.text(dataFile, function (err, text) {
      var columns, data;

      if (err) reject(new Error(err));

      // Replace first line with column names
      columns = ['name', 'count_score', 'gene_length', 'rpkm'];
      text = text.replace(/[^\n]+/, columns.join('\t'));

      data = d3.tsv.parse(text).map(function (row) {
        row.count_score = +row.count_score;
        row.gene_length = +row.gene_length;
        row.rpkm = +row.rpkm;
        return row;
      });

      resolve({
        filename: dataFile,
        data: data
      });
    });
  });
}


/* Drawing */
function Graph(id, height, width) {
  height = height || consts.PLOT_HEIGHT;
  width = width || consts.PLOT_WIDTH;

  this.svg = d3.select(id)
    .append('svg')
    .attr('width', width)
    .attr('height', height)
    .style('position', 'absolute');


  var table = document.getElementById('the-table');
  this.tableSort = new Tablesort(table, { descending: true });

  this.initData();
  this.drawDrawingPad();
  this.initListeners();
}

Graph.prototype = {
  // Create the SVG element that will contain the graph
  initData: function () {
    this.organism = new Organism();
    this.filters = {};
  },

  initListeners: function () {
    d3.select('#scale').on('change', this.draw.bind(this));
  },

  addCell: function (name, cell) {
    var that = this;

    this.organism.addCell(name, cell);
    this.draw();
  },

  removeCell: function (name) {
    this.organism.removeCell(name);
    this.draw();
  },

  clearCells: function () {
    this.organism.reset();
    this.draw();
  },

  addFilter: function (name, filter) {
    this.filters[name] = filter;
    this.draw();
  },

  removeFilter: function (filterId) {
    if (filterId in this.filters) {
      delete this.filters[filterId];
      this.draw();
    }
  },

  clearFilters: function () {
    this.filters = {};
    this.draw();
  },

  passesFilters: function (row, geneName) {
    var ok = true;
    for (var key in this.filters) {
      ok = this.filters[key](row, geneName);
      if (!ok) break;
    }
    return ok;
  },

  getData: function (nonZero) {
    var data = this.organism.asMatrix()
      , allGenes = []
      , filteredGenes = []
      , max = { all: 0, filtered: 0}
      , min = { all: Infinity, filtered: Infinity }

    data.matrix.forEach(function (arr, i) {
      var localMax, localMin;

      if (!(d3.sum(arr))) return;

      allGenes.push({ name: this.organism.genome[i], rpkms: arr });

      localMax = d3.max(arr);
      localMin = d3.min(nonZero ? arr.filter(parseInt) : arr);

      if (localMax > max.all) max.all = localMax;
      if (localMin < min.all) min.all = localMin;

      if (this.passesFilters(arr, this.organism.genome[i])) {
        if (localMax > max.filtered) max.filtered = localMax;
        if (localMin < min.filtered) min.filtered = localMin;
        filteredGenes.push({ name: this.organism.genome[i], rpkms: arr });
      }
    }, this)

    return { cells: data.cells, allGenes: allGenes, filteredGenes: filteredGenes }
  },

  getLines: function (scale, geneMatrix) {
    var NUM_BINS = document.querySelector('#clustering').value

    var all = Array.prototype.concat.apply([], geneMatrix);

    var bins = d3.scale.linear()
      .domain(scale.range())
      .ticks(NUM_BINS)
      .reverse()
      .map(function (px) { return scale.invert(px) })

    var numCells = geneMatrix[0].length;
     
    var genesByBin = bins.reduce(function (acc, bin, binIdx) {acc[binIdx] = []; return acc}, {});
    var binsByGene = geneMatrix.reduce(function (acc, gene, geneIdx) {acc[geneIdx] = []; return acc}, { length: geneMatrix.length });

    for (var cellIdx = 0; cellIdx < numCells; cellIdx++) {
      var cell = geneMatrix
        .map(function (gene, geneIdx) { return { geneIdx: geneIdx, rpkm: gene[cellIdx] } })
        .sort(function (a, b) { return a.rpkm - b.rpkm });
      var binIdx = 0;
      var binLength = bins.length;
      cell.forEach(function (gene) {
        while (true) {
          if (gene.rpkm <= bins[binIdx + 1]) {
            genesByBin[binIdx].push(gene.geneIdx);
            binsByGene[gene.geneIdx].push(binIdx);
            break;
          }
          binIdx += 1;
          if (binIdx >= binLength) throw new Error("Couldn't place all elements in bins.");
        }
      });
    }

    var ret = Array.prototype.reduce.call(binsByGene, function (acc, gene, idx) {
      var seq;
      for (var i = 0; i < numCells - 1; i++) {
        seq = gene.slice(i, i + 2);
        if (!acc[i].hasOwnProperty(seq)) {
          acc[i][seq] = [];
        }
        acc[i][seq].push(idx);
      }
      return acc;
    }, Array.apply(null, { length: numCells - 1 }).map(function () { return {} }));

    var color = d3.scale.linear()
      .range(["blue","red"])
      .domain(ret.reduce(function (acc, seq) {
        var keys = Object.keys(seq);
        if (seq.hasOwnProperty('0,0')) {
          keys.splice(keys.indexOf('0,0'), 1);
        }
        var sorted = keys.map(function (key) {
          return seq[key].length
        }).sort(function (a, b) { return a - b });

        var min = parseInt(sorted.shift(), 10);
        var max = parseInt(sorted.pop(), 10);

        if (min < acc[0]) acc[0] = min;
        if (max > acc[1]) acc[1] = max;

        return acc

      }, [Infinity, 0]));

    return ret.map(function (seqs, idx) {
      return Object.keys(seqs).map(function (key) {
        var points = key.split(',');
        return {
          x1: idx,
          y1: bins[points[0]],
          x2: idx + 1,
          y2: bins[points[1]],
          color: color(seqs[key].length),
          genes: seqs[key]
        }
      });
    });
  },

  draw: function () {
    var log = document.querySelector('#scale_log').checked;

    var data = this.getData(log);

    this.drawTable(data);

    var maxAllRPKM = d3.max(data.allGenes.map(function (gene) { return d3.max(gene.rpkms) }))
      , minAllRPKM = d3.min(data.allGenes.map(function (gene) { return d3.min(log ? gene.rpkms.filter(parseInt) : gene.rpkms) }))
      , maxFilteredRPKM = d3.max(data.filteredGenes.map(function (gene) { return d3.max(gene.rpkms) }))
      , minFilteredRPKM = d3.min(data.filteredGenes.map(function (gene) { return d3.min(log ? gene.rpkms.filter(parseInt) : gene.rpkms) }))

    // Scales
    var y = (log ? d3.scale.log() : d3.scale.linear())
      .domain([maxFilteredRPKM, minFilteredRPKM])
      .range([0 + consts.CELL_PADDING, consts.CELL_HEIGHT - consts.CELL_PADDING])
      .nice()

    if (log) {
      y = y.clamp(true);
    }

    this.y = y;

    var y2 = y.copy().domain([maxAllRPKM, minAllRPKM]).nice();

    var yAxis = d3.svg.axis()
      .scale(y)
      .orient('right')

    var yAxis2 = d3.svg.axis()
      .scale(y2)
      .orient('right')
    
    var tickValues = y.ticks(7);

    var x = function (i) { return 150 + (consts.CELL_WIDTH * i) }
    this.x = x;

    var outlines = {
      'top': y(d3.max(tickValues)),
      'bottom': y(d3.min(tickValues))
    }

    // Draw axes
    this.svg.selectAll('*').remove();

    // Top and bottom line
    this.svg.append('line')
      .attr('x1', x(0))
      .attr('x2', x(data.cells.length - 1))
      .attr('y1', outlines.bottom)
      .attr('y2', outlines.bottom)
      .attr('stroke', 'black')
      .attr('stroke-width', 1)

    this.svg.append('line')
      .attr('x1', x(0))
      .attr('x2', x(data.cells.length - 1))
      .attr('y1', outlines.top)
      .attr('y2', outlines.top)
      .attr('stroke', 'black')
      .attr('stroke-width', 1)

    // Guidelines
    if (false) {
      this.svg.insert('g', ':first-child').selectAll('.guidelines')
        .data(tickValues.slice(1, -1).map(function (yCoord) {
          // x1, y1, x2, y2
          return [x(0), y(yCoord), x(data.cells.length -1), y(yCoord)];
        }))
          .enter()
        .append('line')
        .attr('x1', function (d) { return d[0] })
        .attr('y1', function (d) { return d[1] })
        .attr('x2', function (d) { return d[2] })
        .attr('y2', function (d) { return d[3] })
        .attr('stroke', '#ccc')
    }

    var lineFn = d3.svg.line()
      .x(function (d, i) { return x(i) })
      .y(function (d, i) { return y(d) })
      .interpolate('cardinal')
      .tension(0.85)

    var lines = this.getLines(y, data.filteredGenes.map(function (g) { return g.rpkms }));

    this.svg.selectAll('.gene-group').data(lines)
        .enter()
      .append('g').attr('class', 'gene-group').selectAll('line')
        .data(function (d) { return d.sort(function(a, b) { return a.genes.length - b.genes.length }) })
        .enter()
      .append('line')
      .attr('x1', function (d, i) { return x(d.x1) })
      .attr('y1', function (d, i) { return y(d.y1) })
      .attr('x2', function (d, i) { return x(d.x2) })
      .attr('y2', function (d, i) { return y(d.y2) })
      .attr('stroke', function (d, i) { return d.color })
      .attr('stroke-width', '1')
      .attr('opacity', '.2')

    /*
    this.svg.append('g').selectAll('.gene-path').remove().data(data.filteredGenes)
        .enter()
      .append('path').datum(function (d) { return d.rpkms })
      .classed('gene-path', true)
      .attr('d', lineFn)
      .attr('stroke', 'blue')
      .attr('stroke-width', '1')
      .style('opacity', '.1')
      .attr('fill', 'none')
    */

    this.svg.selectAll('.axis').data(data.cells.map(function (cell) { return { name: cell } }))
        .enter()
      .append('g')
      .attr('class', 'axis')
      .attr('transform', function (d, i) { return 'translate(' + x(i) + ',0)' })
      .each(function (d, i) {
        var first = i === 0
          , last = i === (data.cells.length - 1)

        var axis = d3.select(this)
          .call(
            yAxis
              .orient(first ? 'left': 'right')
              .innerTickSize( (first || last) ? 6 : 12)
          )

        if (! (first || last) ) {
          axis.selectAll('.tick line')
            .attr('transform', 'translate(-6,0)')
          axis.selectAll('.tick text').remove();
        }

        if (log) {
          axis.selectAll('.tick text')
            .text(null)
          .filter(function (d) {
            return d / Math.pow(10, Math.ceil(Math.log(d) / Math.LN10 - 1e-12)) === 1;
          })
          .text(10)
          .append('tspan')
            .attr('dy', '-.7em')
            .text(function (d) { return Math.round(Math.log(d) / Math.LN10) })
        }

        axis
          .append('text')
          .text(function (d) { return d.name })
          .attr('y', outlines.bottom)
          .attr('dy', '.8em')
          .attr('text-anchor', 'start')
          .attr('transform', 'rotate(25, 0, ' + outlines.bottom + ')')

      });

    var brush = d3.svg.brush().y(y2)
      //.on('brush', brushed)

    this.svg.append('g')
      .attr('class', 'brush-axis')
      .call(yAxis2)

    this.svg.append('g')
      .attr('class', 'brush')
      .call(brush)
      .selectAll('rect')
        .attr('width', 50)

    var histogramDataAll = d3.layout.histogram()
      .bins(y2.ticks(50))
      (Array.prototype.concat.apply([], data.allGenes.map(function (gene) { return gene.rpkms })))

    var histogramDataFiltered = d3.layout.histogram()
      .bins(y2.ticks(20))
      (Array.prototype.concat.apply([], data.filteredGenes.map(function (gene) { return gene.rpkms })))

    var histogramX = d3.scale.linear()
      .domain([0, d3.max(histogramDataAll, function (d) { return d.y })])
      .range([0, 50])

    var bar = this.svg.selectAll('.bar')
        .data(histogramDataAll)
      .enter().append('g')
        .attr('class', 'bar')
        .attr('transform', function (d) { return 'translate(' + histogramX(d.y) + ',' + (y2(d.x) - 14) + ')' })

    bar.append('rect')
      .attr('y', 10)
      .attr('width', 10)
      .attr('height', 3)

    var bins = d3.scale.linear().domain(y.range()).ticks(100).reverse().map(function (d) { return y.invert(d) });
    var hist = d3.layout.histogram().bins(bins);

    data.cells.forEach(function (cell, i) {
      var these = data.filteredGenes.map(function (d) { return d.rpkms[i] }).sort();
      var those = data.filteredGenes.map(function (d) { return d.rpkms[i + 1] }).sort();
      var zzz = hist(these);
    }, this);

  },
  
  drawTable: function (data) {
    var that = this;

    var wormbase_url = 'http://www.wormbase.org/search/gene/';

    var toDraw = data.filteredGenes.map(function (gene) {
      var foldChanges = Object.keys(that._foldChanges).map(function (key) {
        var cells = that._foldChanges[key]
          , higher = cells[0]
          , lower = cells[1]

        var orgCells = Object.keys(that.organism.cells)
          , higherIdx = orgCells.indexOf(higher)
          , lowerIdx = orgCells.indexOf(lower)

        if (higherIdx === -1 || lowerIdx === -1) return null;

        var higherVal = gene.rpkms[higherIdx]
          , lowerVal = gene.rpkms[lowerIdx]

        return {
          label: higher + ' -> ' + lower,
          value: higherVal / (lowerVal || .0001)
        }
      }).filter(function (zzz) { return zzz });

      return {
        label: gene.name,
        rpkms: gene.rpkms.slice(),
        avg: d3.mean(gene.rpkms).toFixed(4),
        foldChanges: foldChanges
      }
    });

    d3.select('#table-label').select('div').remove();
    d3.select('#table-label').append('div')
      .html('' + toDraw.length + ' genes ')
      .append('a')
      .html(' (download)')
      .attr('href', '#')
      .on('click', function () {
        var tableHTML = '<table>' + document.getElementById('the-table').innerHTML + '</table>';

        tableHTML = '<html><body>' + '<h1>Results from ' + (new Date()) + '</h1>' + tableHTML + '</body></html>';

        var blob = new Blob([tableHTML], {type: 'text/html'});
        saveAs(blob, 'gene_expression_results.html')
        return false;
      });


    var foldChanges = toDraw[0].foldChanges.map(function (thing) { return thing.label });

    d3.selectAll('.fold-change-header').remove();

    foldChanges.forEach(function (label) {
      d3.select('#table-headers').append('th')
        .attr('class', 'fold-change-header')
        .html(label);
    })

    var select = d3.select('#data-table').selectAll('tr')
      .data(toDraw, function (gene) { return gene.label });

    var lineFn = d3.svg.line()
      .x(function (d, i) { return that.x(i) })
      .y(function (d, i) { return that.y(d) })
      .interpolate('cardinal')
      .tension(0.85)

    var line = null;

    select.enter()
      .append('tr')
      .attr('data-cell', function (d) { return d.label })
      .html(function (d) {
        var url = wormbase_url + d.label;
        var anchor = '<a title="Wormbase" target="_blank" href="' + url + '"><i class="link-icon"></i></a>';

        return '<td>' + anchor + d.label + '</td><td class="mean">' + d.avg + '</td>';
      })
      .on('mouseover', function () {
        var rpkms = d3.select(this).datum().rpkms;
        line = that.svg.append('path').datum(rpkms)
          .classed('gene-path', true)
          .attr('d', lineFn)
          .attr('stroke', 'orange')
          .attr('stroke-width', '2')
          .attr('fill', 'none')

      })
      .on('mouseout', function () {
        line.remove();
      })

    select.select('.mean').html(function (d) { return d.avg });

    select.select('.fold-change-value').remove();
    select.selectAll('.fold-change-value').data(function (d, i) { return d.foldChanges })
      .enter()
        .append('td')
        .attr('class', 'fold-change-value')
        .html(function (d) { return d.value });

    select.exit().remove();

    this.tableSort = new Tablesort(document.getElementById('the-table'), { descending: true });
    //this.tableSort.refresh();
  },
  drawDrawingPad: function () {
    var that = this;

    this._paramFilters = [];
    this._foldChanges = {};

    var el = document.querySelector('#drawing-pad');
    var padSvg = makeDrawingPad(el, DATA_FILES);

    padSvg.on('paramsChange', function (e) {
      var paramsObj = d3.event.detail
      var params;

      for (var param in paramsObj) {

        if (that._paramFilters.indexOf(param) !== -1) continue;

        params = paramsObj[param];

        var higherCell = params.lower;
        var lowerCell = params.higher;
        var foldChange = params.foldChange;

        that._foldChanges[param] = [higherCell, lowerCell];

        that.addFilter(param, function (row) {
          var cells = Object.keys(that.organism.cells);

          var higher = cells.indexOf(higherCell);
          var lower = cells.indexOf(lowerCell);

          if (higher === -1 || lower === -1) {
            return true;
          }

          return (row[higher] / (row[lower] || .0001)) >= foldChange;
        });

        that._paramFilters.push(param);
      }

      that._paramFilters.forEach(function (param) {
        if (!paramsObj.hasOwnProperty(param)) {
          delete that._foldChanges[param];
          that.removeFilter(param);
        }
      });


    });


  }
}


/* Running the thing */

var graph = new Graph('#vis_container');

var li = d3.select('#present-cells').selectAll('li').data(DATA_FILES)
    .enter()
  .append('li')
  .text(function (d) { return d })

li.append('button').text('Add').attr('data-graphed', 'no')
  .on('click', function (d) {
    var that = this;
    if (this.dataset.graphed === 'yes') {
      d3.select(this).text('Add').attr('data-graphed', 'no');
      graph.removeCell(d);
    } else {
      parseFile(d).then(function (dataFile) {
        d3.select(that).text('Remove').attr('data-graphed', 'yes');
        graph.addCell(d, dataFile.data);
      });
    }
  })

/* filters */
var filterCounter = 0;
$('#js-add-filter').on('click', function () {
  var name = 'filter' + filterCounter;
  var $filter = $('#basic-filter')
    .clone()
    .show()
    .prop('id', name)
    .appendTo('#filters');

  filterCounter += 1;

  $filter
    .on('click', '[data-action="save"]', function (e) {
      var quantity = $filter.find('select[name="quantity"]').val()
        , limit = $filter.find('select[name="limit"]').val()
        , value = $filter.find('input').val()

      if (!(quantity && limit && value)) return;

      var fn = function (row) {
        function cmp (d) { return limit === 'gt' ? d > value : d < value }
        if (quantity === 'any') {
          return row.some(cmp);
        } else {
          return row.every(cmp);
        }
      }

      graph.addFilter(name, fn);

      $filter.html('A' + quantity.slice(1) + ' cells have an RPKM ' +
        (limit === 'gt' ? 'greater than ' : 'less than ') +
        value + '. <button data-action="cancel">Remove</button>');
        
    })
    .on('click', '[data-action="cancel"]', function (e) {
      graph.removeFilter(name);
      $filter.remove();
    });
});

