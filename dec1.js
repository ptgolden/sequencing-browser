(function () { 

  var width = 800;
  var height = 200;
  var legendWidth = 100;
  var margin_x = 50;
  var margin_y = 25;
  var axisBuffer = 10;
  var clickable = false;
  var zone = null;
  var newParam = {};
  var finalParams = {};

  var iii = 0;

  var y = d3.scale.log()
      //	.base(100)
        .domain([1,10000])
        .range([0, height-margin_y]);

  function dispatchEvent(svgContainer) {
    var paramEvent = new CustomEvent('paramsChange', { detail: finalParams });
    svgContainer[0][0].dispatchEvent(paramEvent);
  }

  //Set up the canvas on which axes will be draw 
  //  and parameters can be specified by user
  function makeDrawingPad(container, samples) {
    var paramLine, paramGroup;

    sampleNumb = samples.length;
    sensitiveZone = {};
    axisPosition = {};
    
    var svgContainer = d3.select(container).append("svg")
      .attr("width", width+2*margin_x)
      .attr("height", height+2*margin_y)
      .attr("class", "padbg");

    /*
    var cellguys = svgContainer.selectAll('.cellguy').data(samples, function (d) { return d });
    */
    
    //Draw vertical lines and add titles for each sample
    for (var i = 0; i<sampleNumb; i++) {
      svgContainer.append('g');

      axisPosition[samples[i]] = margin_x + i*(width/(sampleNumb-1));
      svgContainer.append("line")
        .attr("x1", axisPosition[samples[i]])
        .attr("y1", 2*margin_y)
        .attr("x2", axisPosition[samples[i]])
        .attr("y2", height+margin_y)
        .attr("class", "axisLine")
        .attr("stroke", "darkgray");
        
      svgContainer.append("text")
        .attr("text-anchor", "middle")
        .attr("x", axisPosition[samples[i]])
        .attr("y", margin_y*1.25)
        .attr("class", "notSelectable")
        .text(samples[i]);	

      //Specify the regions that should register as axes when clicked.
      svgContainer.append("rect")
        .attr("x", axisPosition[samples[i]]-axisBuffer)
        .attr("y", 2*margin_y-axisBuffer)
        .attr("width", 2*axisBuffer)
        .attr("height", height-margin_y+2*axisBuffer)
        .style("fill-opacity", 0)
        .attr("id", samples[i])
        .attr("class", "sensitiveZone")
        .style("cursor", "crosshair")

        //When moused over, mouse registers as being in a clickable region,
        //   and the name of the sample is logged temporarily
        .on("mouseover", function(){
          clickable=true;
          zone = this.id;
        })

        //When moused out, sample and clickability are reset to null.
        .on("mouseout", function(){
          clickable=false;
          zone=null;
        });
          
    
      //Make a vector of parameters (empty, so far)
      paramCoords = [];
    
      var drawParam = svgContainer
        .on("mousedown", mousedown)
        .on("mouseup", mouseup);
        
      var mousedown = function() {
        var m = d3.mouse(this);
        if (clickable) {

          iii += 1;

          paramGroup = svgContainer.append('g')
            .attr('id', 'param-group' + iii);

          paramLine = paramGroup.append("line")
            .attr("x1", axisPosition[zone])
            .attr("y1", m[1])
            .attr("x2", m[0])
            .attr("y2", m[1])
            .attr("class", "paramLine");
          newParam["x1"]= zone;
          newParam["y1"]=m[1];
          //drop a circle or figure out a marker start situation
        
        //As mouse is moving, calculate the foldchange for current line
        var m = d3.mouse(this);
        if (newParam["y1"] < m[1]) {
          var FC = d3.round(y.invert(m[1]-newParam["y1"]),2);
        } else {
          var FC = d3.round(y.invert(newParam["y1"]-m[1]),2);
        }
        
        FoChLine = paramGroup.append("text")
          .attr("text-anchor", "middle")
          .attr("class", "FoCh")
          .attr("fill", "red")
          .attr("x", m[0])
          .attr("y", m[1]-margin_y)
          .attr("class", "notSelectable")			
          .text(FC);
          
          svgContainer.on("mousemove", mousemove);			
        };

        paramGroup.on('dblclick', function () {
          var grp = d3.select(this);
          var id = grp.attr('id');

          delete finalParams[id];

          grp.remove();
          dispatchEvent(svgContainer);
        });

      }

      var mousemove = function() {
        //As mouse is moving, calculate the foldchange for current line
        var m = d3.mouse(this);
        if (newParam["y1"] < m[1]) {
          var FC = d3.round(y.invert(m[1]-newParam["y1"]),2);
        } else {
          var FC = d3.round(y.invert(newParam["y1"]-m[1]),2);
        }

        //Shift endpoint of line over so it doesnt cover up the axes
        if (m[0]< axisPosition[newParam["x1"]]) {
          spacer = 3;
        } else {spacer = -3};
        
        //Draw line
        paramLine.attr("x2", m[0]+spacer)
          .attr("y2", m[1]);
        
        //Append fold change annotation
        FoChLine.attr("x", (axisPosition[newParam["x1"]]+m[0])/2)
          .attr("y", ((newParam["y1"]+m[1])/2)-margin_y)			
          .text(FC);
      }
      
      var mouseup = function(){

        var paramEvent;
        var id;

        if (clickable && zone != newParam["x1"]) {

          id = paramLine[0][0].parentNode.getAttribute('id');
          
          //remove the newParam line, and make a permanent line
          paramLine.attr("x2", axisPosition[zone])
            .attr("y2", d3.mouse(this)[1]);
          newParam["x2"] = zone;
          newParam["y2"] = d3.mouse(this)[1];
          paramCoords.push(newParam);
          var wordyParam = {};
          if (newParam["y1"] > newParam["y2"]) {
            wordyParam["higher"] = newParam["x1"];
            wordyParam["lower"] = newParam["x2"];
            wordyParam["foldChange"] = d3.round(y.invert(newParam["y1"]-newParam["y2"]),2);
          } else {
            wordyParam["higher"] = newParam["x2"];
            wordyParam["lower"] = newParam["x1"];
            wordyParam["foldChange"] = d3.round(y.invert(newParam["y2"]-newParam["y1"]),2);
          }
          finalParams[id] = wordyParam;
          newParam = {};
          wordyParam = {};
          svgContainer.on("mousemove", null);

          dispatchEvent(svgContainer);
        }
        else {
          //Only erase line if a new one has been started
          if (newParam["x1"]) {
            paramGroup.remove();
            svgContainer.on("mousemove", null);
          }
        }
      }		

    }

    return svgContainer;
  };

  window.makeDrawingPad = makeDrawingPad;
})();
