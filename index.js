// Constants
var consts = {
  // Height of total plot area
  PLOT_WIDTH: 1200,
  PLOT_HEIGHT: 800,

  // Height and width of the cell graphs
  CELL_HEIGHT: 500,
  CELL_WIDTH: 300,

  // Padding above and below cell graphs
  CELL_PADDING: 30
}


//
var DATA_FILES = [
  'data/193_4cell_RPKM.txt',
  'data/194_4cell_RPKM.txt',
  'data/195_4cell_RPKM.txt',
  'data/196_4cell_RPKM.txt'
];

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

function drawGraph() {
  var svg = d3.select('#vis_container')
    .append('svg')
    .attr('width', consts.PLOT_WIDTH)
    .attr('height', consts.PLOT_HEIGHT)
}

drawGraph();
Promise.all([
  parseFile(DATA_FILES[0]),
  parseFile(DATA_FILES[1]),
  parseFile(DATA_FILES[2]),
  parseFile(DATA_FILES[3])
]).then(drawSamples);

var filters = [];

filters.push(function (gene) {
  return gene.rpkms.some(function (d) { return d > 1500 });
});

function drawSamples(parsedSamples) {
  var genePresences = parsedSamples.reduce(function (acc, sample) {
    sample.data.forEach(function (gene) {
      if (!acc.hasOwnProperty(gene.name)) acc[gene.name] = [];
      acc[gene.name].push(gene.rpkm)
    });
    return acc;
  }, {});

  var presentGenes = Object.keys(genePresences).filter(function (geneName) {
    return d3.sum(genePresences[geneName]);
  });
  
  var geneLines = presentGenes.map(function (geneName) {
    return {
      name: geneName,
      rpkms: genePresences[geneName]
    }
  });

  filters.forEach(function (filter) {
    geneLines = geneLines.filter(filter);
  });

  var maxRpkm = d3.max(geneLines, function (gene) { return d3.max(gene.rpkms) });
  var minRpkm = d3.min(geneLines, function (gene) { return d3.min(gene.rpkms) });

  var y = d3.scale.linear()
    .domain([maxRpkm, minRpkm])
    .range([0 + consts.CELL_PADDING, consts.CELL_HEIGHT - consts.CELL_PADDING])
    .nice()

  var tickValues = y.ticks(7).concat(y.domain());
  var yAxis = d3.svg.axis()
    .scale(y)
    .orient('right')
    .tickValues(tickValues)

  var x = function (i) { return 50 + (300 * i) }

  var line = d3.svg.line()
    .x(function (d, i) { return x(i) })
    .y(function (d, i) { return y(d) })
    .interpolate('cardinal')
    .tension(0.85)

  d3.select('svg').selectAll('path').data(geneLines)
      .enter()
    .append('path').datum(function (d) { return d.rpkms })
    .attr('d', line)
    .attr('stroke', 'blue')
    .attr('stroke-width', '1')
    .style('opacity', '.4')
    .attr('fill', 'none')

  var outlines = {
    'top': y(d3.max(tickValues)),
    'bottom': y(d3.min(tickValues))
  }

  d3.select('svg').selectAll('.axis').data(parsedSamples.map(function (d) { return d.filename }))
      .enter()
    .append('g')
    .attr('class', 'axis')
    .attr('transform', function (d, i) { return 'translate(' + x(i) + ',0)' })
    .each(function (d, i) {
      var first = i === 0
        , last = i === (parsedSamples.length - 1)

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

      axis
        .append('text')
        .text(function (d) { return d })
        .attr('y', outlines.bottom)
        .attr('dy', '.8em')
        .attr('text-anchor', 'start')
        .attr('transform', 'rotate(25, 0, ' + outlines.bottom + ')')

    })

  d3.select('svg').insert('g', ':first-child').selectAll('.guidelines')
    .data(tickValues.map(function (yCoord) {
      // x1, y1, x2, y2
      return [x(0), y(yCoord), x(parsedSamples.length -1), y(yCoord)];
    }))
      .enter()
    .append('line')
    .attr('x1', function (d) { return d[0] })
    .attr('y1', function (d) { return d[1] })
    .attr('x2', function (d) { return d[2] })
    .attr('y2', function (d) { return d[3] })
    .attr('stroke', '#ccc')

  d3.select('svg').append('line')
    .attr('x1', x(0))
    .attr('x2', x(parsedSamples.length - 1))
    .attr('y1', outlines.bottom)
    .attr('y2', outlines.bottom)
    .attr('stroke', 'black')
    .attr('stroke-width', 2)

  d3.select('svg').append('line')
    .attr('x1', x(0))
    .attr('x2', x(parsedSamples.length - 1))
    .attr('y1', outlines.top)
    .attr('y2', outlines.top)
    .attr('stroke', 'black')
    .attr('stroke-width', 2)
}
