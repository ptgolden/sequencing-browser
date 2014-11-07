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


/* Data structures */
function Organism() {
  this.reset();
}

Organism.prototype = {
  reset: function () {
    this.cells = {};
    this.genome = null;
  },

  addCell: function (name, data) {
    var genomeCheck
      , rpkms = []

    if (!this.genome) {
      this.genome = data.map(function (gene) { return gene.name });
    }

    genomeCheck = this.genome.slice();

    data.forEach(function (gene) {
      var check = genomeCheck.shift();
      if (gene.name !== check) throw Error('Cell\'s genome not identical to organism\'s genome');
      rpkms.push(gene.rpkm);
    });

    this.cells[name] = rpkms;
  },

  removeCell: function (name) {
    delete this.cells[name];
  },

  getGeneRPKMS: function (geneName) {
    var idx = this.genome.indexOf(geneName)
      , genes = {}

    if (idx === -1) throw Error('' + geneName + ' is not in this organism\'s genome.');

    for (var cell in this.cells) {
      genes[cell] = this.cells[cell][idx];
    }

    return genes;
  },

  asMatrix: function () {
    var cellNames = Object.keys(this.cells)
      , matrix

    matrix = this.genome.map(function (gene, i) {
      return cellNames.map(function (cellName) {
        return this.cells[cellName][i];
      }, this);
    }, this);

    return { cells: cellNames, matrix: matrix }
  }
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

    var container = document.querySelector('#present-cells');

    var el = document.createElement('li');
    el.classList.add('present-cell');
    el.textContent = name;

    var button = document.createElement('button');
    button.textContent = 'Remove';
    button.addEventListener('click', function () {
      that.removeCell(name);
      container.removeChild(el);
    }, false);
    el.appendChild(button);
    container.appendChild(el);

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

  addFilter: function (filter) {
    this.filters[filter.id] = filter;
    this.draw();
  },

  removeFilter: function (filterId) {
    delete this.filters[filterId];
    this.draw();
  },

  clearFilters: function () {
    this.filters = {};
    this.draw();
  },

  get data() {
    var data = this.organism.asMatrix()
      , genes = []

    // Filter out all non-zero numbers
    data.matrix = data.matrix.filter(function (arr) { return d3.sum(arr) });

    data.matrix.forEach(function (arr, i) {
      if (!(d3.sum(arr))) return;
      genes.push({ name: this.organism.genome[i], rpkms: arr });
    }, this)

    // Apply filters
    // TODO

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
    if (!log) {
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

    var min = tickValues[0];

    this.svg.selectAll('.gene-path').remove().data(data.genes)
        .enter()
      .append('path').datum(function (d) {
        if (!log) return d.rpkms;

        return d.rpkms.map(function (rpkm) { return rpkm || min });
      })
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

parseFile(DATA_FILES[0]).then(function (dataFile) {
  graph.addCell(dataFile.filename, dataFile.data);
});
parseFile(DATA_FILES[1]).then(function (dataFile) {
  graph.addCell(dataFile.filename, dataFile.data);
});
parseFile(DATA_FILES[2]).then(function (dataFile) {
  graph.addCell(dataFile.filename, dataFile.data);
});
