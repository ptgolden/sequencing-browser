// jshint laxcomma: true, asi: true

/* Constants */

var consts = {
  // Height of total plot area
  PLOT_WIDTH: 1200,
  PLOT_HEIGHT: 800,

  // Height and width of the cell graphs
  CELL_HEIGHT: 500,
  CELL_WIDTH: 300,

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
    .attr('height', height);

  this.initData();
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

  passesFilters: function (row) {
    var ok = true;
    for (var key in this.filters) {
      ok = this.filters[key](row);
      if (!ok) break;
    }
    return ok;
  },

  get data() {
    var data = this.organism.asMatrix()
      , genes = []

    data.matrix.forEach(function (arr, i) {
      if (!(d3.sum(arr))) return;
      if (this.passesFilters(arr)) {
        genes.push({ name: this.organism.genome[i], rpkms: arr });
      }
    }, this)

    return { cells: data.cells, genes: genes }
  },

  draw: function () {
    var data = this.data;

    var log = document.querySelector('#scale_log').checked;

    var maxRPKM = d3.max(data.genes.map(function (g) { return d3.max(g.rpkms) }))
      , minRPKM = d3.min(data.genes.map(function (g) {
        if (log) {
          return d3.min(g.rpkms.filter(function (d) { return d }));
        }
        return d3.min(g.rpkms)
      }))

    // Scales
    var y = (log ? d3.scale.log() : d3.scale.linear())
      .domain([maxRPKM, minRPKM])
      .range([0 + consts.CELL_PADDING, consts.CELL_HEIGHT - consts.CELL_PADDING])
      .nice()

    if (log) {
      y = y.clamp(true);
    }

    var yAxis = d3.svg.axis()
      .scale(y)
      .orient('right')
      //.tickValues(tickValues)
    
    var tickValues = y.ticks(7);

    var x = function (i) { return 50 + (300 * i) }

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

    var lineFn = d3.svg.line()
      .x(function (d, i) { return x(i) })
      .y(function (d, i) { return y(d) })
      .interpolate('cardinal')
      .tension(0.85)

    this.svg.selectAll('.gene-path').remove().data(data.genes)
        .enter()
      .append('path').datum(function (d) { return d.rpkms })
      .classed('gene-path', true)
      .attr('d', lineFn)
      .attr('stroke', 'blue')
      .attr('stroke-width', '1')
      .style('opacity', '.1')
      .attr('fill', 'none')
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

