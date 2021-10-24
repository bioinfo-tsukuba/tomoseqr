"use strict";

Shiny.addCustomMessageHandler("initData", InitData);
Shiny.addCustomMessageHandler("existingData", ExistingData);

const mainCanvas = document.getElementById("sample");
const mainCtx = mainCanvas.getContext('2d');

class MaskFigure {
  constructor(initData, maskDim, xRatio, yRatio, alongAxis) {
    this.position = 0;
    this.dataMatrix = initData;
    this.xRatio = xRatio;
    this.yRatio = yRatio;
    this.xLen = maskDim[0] * this.xRatio;
    this.yLen = maskDim[1] * this.yRatio;
    this.numberOfPages = maskDim[2];
    this.isDrawing = false;
    this.figureSet = [];
    this.alongAxis = alongAxis;
    this.mousedown = (e) => {this._MouseDown(e)};
    this.mousemove = (e) => {this._MouseMove(e)};
    this.mouseup = (e) => {this._MouseUp(e)};
  }

  SendToR() {
    Shiny.setInputValue("mask_matrix", this.dataMatrix.flat(Infinity));
    Shiny.setInputValue("along_axis", this.alongAxis);
    Shiny.setInputValue("mask_dim", [this.dataMatrix.length, this.dataMatrix[0].length, this.dataMatrix[0][0].length])
  }

  _MouseDown(e) {
    const mX = Math.floor(e.offsetX / this.xRatio);
    const mY = Math.floor((this.yLen - e.offsetY) / this.yRatio);
    if (pencil == true) {
      mainCtx.fillStyle = "black";
      mainCtx.fillRect(mX*this.xRatio, -1 * (mY - (this.yLen / this.yRatio) + 1)*this.yRatio, this.xRatio, this.yRatio);
      this.dataMatrix[mX][mY][this.position] = 1;
    } else {
      mainCtx.fillStyle = "white";
      mainCtx.fillRect(mX*this.xRatio, -1 * (mY - (this.yLen / this.yRatio) + 1)*this.yRatio, this.xRatio, this.yRatio);
      this.dataMatrix[mX][mY][this.position] = 0;
    }
    this.isDrawing = true;
    this.addGrid("white");
    this.addGrid("gray");
  }

  _MouseMove(e) {
    const mX = Math.floor(e.offsetX / this.xRatio);
    const mY = Math.floor((this.yLen - e.offsetY) / this.yRatio);
    document.getElementById("txtY").innerHTML = 'Position: (' + (mX + 1) + ', ' + (mY + 1) + ')';
    if (this.isDrawing == true) {
      if (pencil == true) {
      mainCtx.fillStyle = "black";
      mainCtx.fillRect(mX*this.xRatio, -1 * (mY - (this.yLen / this.yRatio) + 1)*this.yRatio, this.xRatio, this.yRatio);
      this.dataMatrix[mX][mY][this.position] = 1;
    } else {
      mainCtx.fillStyle = "white";
      mainCtx.fillRect(mX*this.xRatio, -1 * (mY - (this.yLen / this.yRatio) + 1)*this.yRatio, this.xRatio, this.yRatio);
      this.dataMatrix[mX][mY][this.position] = 0;
    }
    this.addGrid("white");
    this.addGrid("gray");
    }
  }

  _MouseUp(e) {
    this.isDrawing = false;
    this.addGrid("white");
    this.addGrid("gray");
    this.SendToR();
  }

  RemoveListener() {
    mainCanvas.removeEventListener('mousedown', this.mousedown);
    mainCanvas.removeEventListener('mousemove', this.mousemove);
    mainCanvas.removeEventListener('mouseup', this.mouseup);
  }

  InitCanvas() {
    mainCanvas.setAttribute("width", this.xLen);
    mainCanvas.setAttribute("height", this.yLen);
    mainCtx.fillStyle = "white";
    mainCtx.clearRect(0, 0, 9999, 9999)
    mainCtx.fillRect(0, 0, this.xLen, this.yLen);
    document.getElementById("nowPage").innerHTML = "Page: " + (this.position + 1) + ' / ' + this.numberOfPages;
    document.getElementById("txtY").innerHTML = "Position: ";

    mainCanvas.addEventListener('mousedown', this.mousedown);
    mainCanvas.addEventListener('mousemove', this.mousemove);
    mainCanvas.addEventListener('mouseup', this.mouseup);

  }
  addGrid(color) {
    const x = this.xLen;
    const y = this.yLen;
    mainCtx.strokeStyle = color;
    mainCtx.lineWidth = 1;
    mainCtx.beginPath();
    for (let i = 0; i < x + 1; i += this.xRatio) {
      mainCtx.moveTo(i, 0);
      mainCtx.lineTo(i, y);
    }
    for (let i = 0; i < y + 1; i += this.yRatio) {
      mainCtx.moveTo(0, i);
      mainCtx.lineTo(x, i);
    }
    mainCtx.stroke();
  }

  ChangeFigure(step) {
    this.figureSet[this.position] = mainCtx.getImageData(0, 0, this.xLen, this.yLen);
    if((this.position + step) < 0) {
      alert("This is the first page of figure!");
      return 1;
    }
    if((this.position + step) >= this.numberOfPages) {
      alert("This is the last page of figure!");
      return 1;
    }
    this.position += step;
    mainCtx.putImageData(this.figureSet[this.position], 0, 0);
    document.getElementById("nowPage").innerHTML = "Page: " + (this.position + 1) + ' / ' + this.numberOfPages;
  }


}

function DeNovoFigure(initData, maskDim, xRatio, yRatio, alongAxis) {
  const initCanvas = new MaskFigure(initData, maskDim, xRatio, yRatio, alongAxis);
  initCanvas.InitCanvas();
  initCanvas.addGrid();
  for(let i = 0; i < maskDim[2]; i++) {
    initCanvas.figureSet.push(mainCtx.getImageData(0, 0, initCanvas.xLen, initCanvas.yLen));
  }
  return initCanvas;
}

/*
Reflect an information of existing mask matrix to canvas.
This function will used in making figureSet.
*/
function DrawXY(dataMatrix, z, xRatio, yRatio) {
  const xLen = dataMatrix.length;
  const yLen = dataMatrix[0].length;
  const tmpCanvas = document.getElementById("tmpCanvas");
  const tmpCtx = tmpCanvas.getContext('2d');
  tmpCanvas.setAttribute("width", xLen * xRatio);
  tmpCanvas.setAttribute("height", yLen * yRatio);
  tmpCtx.fillStyle = "white";
  tmpCtx.clearRect(0, 0, 9999, 9999)
  tmpCtx.fillRect(0, 0, xLen * xRatio, yLen * yRatio);
  for(let i = 0; i < xLen; i++) {
    for(let j = 0; j < yLen; j++) {
      if(dataMatrix[i][j][z] == 1) {
        tmpCtx.fillStyle = "black";
        tmpCtx.fillRect(i * xRatio, (yLen - j - 1) * yRatio, xRatio, yRatio);
      }
    }
  }
  return tmpCtx.getImageData(0, 0, xLen * xRatio, yLen * yRatio);
}

function FigureFromData(dataMatrix, xRatio, yRatio, alongAxis){
  const xLen = dataMatrix.length;
  const yLen = dataMatrix[0].length;
  const zLen = dataMatrix[0][0].length;
  const initCanvas = new MaskFigure(dataMatrix, [xLen, yLen, zLen], xRatio, yRatio, alongAxis);
  initCanvas.InitCanvas();
  initCanvas.addGrid();
  for(let z = 0; z < zLen; z++) {
    initCanvas.figureSet.push(DrawXY(dataMatrix, z, xRatio, yRatio));
  }
  initCanvas.addGrid("white");
  initCanvas.addGrid("gray");
  return initCanvas;
}


let initCanvas;
let pencil = true;
function InitData(message) {
  if(typeof initCanvas != 'undefined') {
    initCanvas.RemoveListener();
  }
  initCanvas = DeNovoFigure(message.dataMatrix, message.maskDim, message.ratio[0], message.ratio[1], message.alongAxis);
  initCanvas.SendToR();
}

function ExistingData(message) {
  if(typeof initCanvas != 'undefined') {
    initCanvas.RemoveListener();
  }
  initCanvas = FigureFromData(message.dataMatrix, message.ratio[0], message.ratio[1], message.alongAxis);
  mainCtx.putImageData(initCanvas.figureSet[0], 0, 0);
  initCanvas.addGrid("gray");
  initCanvas.SendToR();
}

function next_click(){
  initCanvas.ChangeFigure(1);
  initCanvas.addGrid("white");
  initCanvas.addGrid("gray");
}
function back_click(){
  initCanvas.ChangeFigure(-1);
  initCanvas.addGrid("white");
  initCanvas.addGrid("gray");
}

function SetPencil(){
  pencil = true;
}
function SetEraser(){
  pencil = false;
}