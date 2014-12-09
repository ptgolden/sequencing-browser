// jshint laxcomma: true, asi: true
(function (global) {
  "use strict";

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
        if (gene.name !== check) throw new Error('Cell\'s genome not identical to organism\'s genome');
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

      if (idx === -1) throw new Error('' + geneName + ' is not in this organism\'s genome.');

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

  global.Organism = Organism;
})(this);
