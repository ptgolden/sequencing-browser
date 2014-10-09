// Constants
var consts = {
  // Height of total plot area
  PLOT_WIDTH: 1800,
  PLOT_HEIGHT: 1000,

  // Height and width of the cell graphs
  CELL_HEIGHT: 500,
  CELL_WIDTH: 300,

  // Padding between cell graphs
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

      resolve(data);
    });
  });
}

function drawGraph() {
  var svg = d3.select('#vis_container')
    .append('svg')
    .attr('width', consts.PLOT_WIDTH)
    .attr('height', consts.PLOT_HEIGHT)
    .style('background', '#f0f0f0')
}

drawGraph();
Promise.all([
  parseFile(DATA_FILES[0]),
  parseFile(DATA_FILES[1]),
  parseFile(DATA_FILES[2]),
  parseFile(DATA_FILES[3])
]).then(drawSamples);

function drawSamples(parsedSamples) {
  var maxRpkm = d3.max(parsedSamples, function (genes) {
    return d3.max(genes, function (gene) {
      return gene.rpkm;
    });
  });

  var genePresences = parsedSamples.reduce(function (acc, sample) {
    sample.forEach(function (gene) {
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
      coords: genePresences[geneName].reduce(function (acc, el, i, all) {
        if (i === 0) return acc;
        acc.push([all[i-1], el]);
        return acc;
      }, [])
    }
  });

  var y = d3.scale.linear()
    .domain([maxRpkm, 0])
    .range([0, consts.PLOT_HEIGHT])

  var cells = d3.select('svg').selectAll('.cells').data(parsedSamples);

  cells.enter().append('g')
    .attr('class', 'cells')
    .attr('transform', function (d, i) {
      return 'translate(' + (300 * (i + 1)) + ',0)';
    });

  cells.selectAll('circle')
    .data(function (d) {
      return d.filter(function (gene) { return presentGenes.indexOf(gene.name) !== -1 });
    })
    .enter()
    .append('circle')
    .attr('r', '1')
    .attr('cx', 0)
    .attr('cy', function (d) { return y(d.rpkm) })

  var lines = d3.select('svg').selectAll('.lines').data(geneLines)
        .enter()
      .append('g').attr('class', 'lines');

  lines.selectAll('line').data(function (d) { return d.coords })
      .enter()
    .append('line')
    .attr('x1', function (d, i) { return 300 * (i + 1) })
    .attr('y1', function (d) { return y(d[0]) })
    .attr('x2', function (d, i) { return 300 * (i + 2) })
    .attr('y2', function (d) { return y(d[1]) })
    .attr('stroke', 'blue')
    .style('opacity', '.1')
}
