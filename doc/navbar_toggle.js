// JavaScript Document
/* 
The Visibility Toggle
Copyright 2003 by Sim D'Hertefelt 
www.interactionarchitect.com 
info@interactionarchitect.com
*/

var hotspots = document.getElementsByName('hotspot');
var toggles = document.getElementsByName('toggle');

function visibilitytoggle()
{
  for (var i = 0; i < hotspots.length; i++)
  {
  hotspots[i].someProperty = i;
  hotspots[i].onclick = function() {toggle(this.someProperty)};
  }

  for (var i = 0; i < toggles.length; i++)
  {
  toggles[i].style.display = 'none';
  }
}

function toggle(i)
{
  if (toggles[i].style.display == 'none')
  {
   hideall()     
   toggles[i].style.display = ''
  }
  else
  toggles[i].style.display = 'none'
} 

function showall()
{
  for (var i = 0; i < toggles.length; i++)
  {
  toggles[i].style.display = '';
  }
}

function hideall()
{
  for (var i = 0; i < toggles.length; i++)
  {
  toggles[i].style.display = 'none';
  }
}

