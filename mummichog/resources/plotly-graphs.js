
//***** MZ scatter plot
function makeMZUserInputPlot() {
    //Plotly.d3.json("result.json", function(data){ processMZData(data.userData_original, cutoff) } );
    processMZData(data.userData_original, data.meta_data.input_parameters.cutoff)
};

function processMZData(allRows, cutoff) {

  console.log(allRows);
  var xTrace1 = [], yTrace1 = []
  var xTrace2 = [], yTrace2 = []

  for (var i=0; i<allRows.length; i++) {
    row = allRows[i];
    if (row['p_value'] < parseFloat(cutoff)){
        xTrace1.push( parseFloat(row['mz']) );
        yTrace1.push( - Math.log10(row['p_value']) );
    }else{
        xTrace2.push( parseFloat(row['mz']) );
        yTrace2.push( - Math.log10(row['p_value']) );
    }
  }

  makeMZPlotly( xTrace1, yTrace1, xTrace2, yTrace2, cutoff);
}

function makeMZPlotly( xTrace1, yTrace1, xTrace2, yTrace2, cutoff){
  var plotDiv = document.getElementById("plot");
  var yLineValue = - Math.log10(cutoff)
  var maxMz = Math.max(...xTrace1, ...xTrace2)
    console.log('maxMz=' + maxMz)
  var trace1 = {
    mode: 'markers',
    x: xTrace1,
    y: yTrace1,
    name: 'Significant, p < ' + cutoff + ', (' + xTrace1.length + ')',
    type: 'scatter'
  };

  var trace2 = {
    mode: 'markers',
    x: xTrace2,
    y: yTrace2,
    name: 'p >= ' + cutoff + ', (' + xTrace2.length + ')',
    type: 'scatter'
  };

  var trace3 = {
    mode: 'lines',
    x: [0, maxMz],
    y: [yLineValue, yLineValue],
    line: {
        dash: 'dot',
        width: 2,
        color: 'rgb(128,0,128)'
    },
    name: 'p'
  };

  var data = [trace1, trace2, trace3];

  var layout = {
      xaxis: {
         autorange: true,
         title: 'm/z'
      },
      yaxis: {
          autorange: true,
          title: '-log10 p-value'
      },
      title:'User Input - m/z'
    };

  Plotly.newPlot('mz_user_input', data, layout);
};


//*** Ret Time scatter plot
function makeRetTimeUserInputPlot() {
  //Plotly.d3.csv("userInputData.csv", function(data){ processRetTimeData(data) } );
    processRetTimeData(data.userData_original, data.meta_data.input_parameters.cutoff)
};

function processRetTimeData(allRows, cutoff) {

  console.log(allRows);
  var xTrace1 = [], yTrace1 = []
  var xTrace2 = [], yTrace2 = []

  for (var i=0; i<allRows.length; i++) {
    row = allRows[i];
    if (row['p_value'] < parseFloat(cutoff)){
        xTrace1.push( row['retention_time'] );
        yTrace1.push( - Math.log10(row['p_value']) );
    }else{
        xTrace2.push( row['retention_time'] );
        yTrace2.push( - Math.log10(row['p_value']) );
    }
  }

  makeRetTimePlotly( xTrace1, yTrace1, xTrace2, yTrace2, cutoff);
}

function makeRetTimePlotly( xTrace1, yTrace1, xTrace2, yTrace2, cutoff){
  var plotDiv = document.getElementById("plot");
  var yLineValue = - Math.log10(cutoff)
  var maxRetTime = Math.max(...xTrace1, ...xTrace2)

  var trace1 = {
    mode: 'markers',
    x: xTrace1,
    y: yTrace1,
    name: 'Significant, p < ' + cutoff + ', (' + xTrace1.length + ')',
    type: 'scatter'
  };

  var trace2 = {
    mode: 'markers',
    x: xTrace2,
    y: yTrace2,
    name: 'p >= ' + cutoff + ', (' + xTrace2.length + ')',
    type: 'scatter'
  };

  var trace3 = {
    mode: 'lines',
    x: [0, maxRetTime],
    y: [yLineValue, yLineValue],
    line: {
        dash: 'dot',
        width: 2,
        color: 'rgb(128,0,128)'
    },
    name: 'p'
  };

  var data = [trace1, trace2, trace3];

  var layout = {
      xaxis: {
         autorange: true,
         title: 'Retention time'
      },
      yaxis: {
          autorange: true,
          title: '-log10 p-value'
      },
      title:'User Input - Retention Time'
    };


  Plotly.newPlot('retention_time_input', data, layout);
};


/// ** Pathway bar plot

function makePathwayPlot() {
    //Plotly.d3.csv("mcg_pathwayanalysis_myResult.csv", function(data){ processPathwayData(data) } );
    processPathwayData(data.result_pathwayAnalysis, data.meta_data.significance_cutoff)
};

function processPathwayData(allRows, cutoff) {

  console.log(allRows);
  var xTrace1 = [], yTrace1 = [], xLineTrace = [];

  for (var i=0; i<allRows.length; i++) {
    row = allRows[i];
    if (row['pathway_p_value'] < parseFloat(cutoff)){
        yTrace1.push( row['name'] );
        xTrace1.push( - Math.log10(row['pathway_p_value']) );
        xLineTrace.push (1.301)
    }
  }

  console.log( 'X tarce 1',xTrace1, 'Y trace 1', yTrace1);
  makePathwayPlotly( xTrace1, yTrace1, xLineTrace);
}

function makePathwayPlotly( xTrace1, yTrace1, xLineTrace){
  var plotDiv = document.getElementById("plot");
  var trace1 = {
    x: xTrace1,
    y: yTrace1,
    type: 'bar',
    orientation: 'h',
    marker: {
        color: 'rgba(50,171,96,0.6)'
    },
    transforms: [{
        type: 'sort',
        target: 'y',
        order: 'descending'
    }]
  };

  var trace2 = {
    x: xLineTrace,
    y: yTrace1,
    mode: 'lines',
    orientation: 'h',
    line: {
        dash: 'dot',
        width: 2,
        color: 'rgb(128,0,128)'
    }
  };

  var data = [trace1, trace2];

  var layout = {
      yaxis: {
         autorange: true,
         title: 'Pathways'
      },
      xaxis: {
          autorange: true,
          title: '-log10 p-value'
      },
      title:'Pathway Analysis',
      margin: {
        l: 320,
        r: 20,
        t: 100,
        b: 70
      },
      showlegend: false
    };


  Plotly.newPlot('pathwayBarPlot', data, layout);
};

