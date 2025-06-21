// JavaScript Document
/*
The Visibility Toggle
Copyright 2003 by Sim D'Hertefelt
www.interactionarchitect.com
info@interactionarchitect.com
*/
// Bits to click
var hotspots = document.getElementsByName('hotspot');
// Menus that appear when hotspot is clicked
var toggles = document.getElementsByName('toggle');

var searchActivatedFlag = false;


function visibilitytoggle()
{
  for (var i = 0; i < hotspots.length; i++)
  {
    hotspots[i].someProperty = i;
    hotspots[i].onclick = function() {toggle(this.someProperty)};
  }
}

function toggle(i)
{
  if (toggles[i].style.display == 'none')
  {
    hideall()

    toggles[i].style.display = ''
    if(i == 3) document.getElementById('search_wrapper_outer').style.display = '';
    if(i != 3) document.getElementById('extended_menu_outer').style.display = ''; // 4 is the search bar, this is pretty hacky, need to reconsider

    hotspots[i].style.backgroundColor = '#e5e5e5';
    hotspots[i].style.borderBottom = '1px solid #e5e5e5';
    hotspots[i].style.borderRight = '1px solid  #d5d5d5';
    hotspots[i].style.borderLeft = '1px solid  #d5d5d5';

    if(i==3 && searchActivatedFlag == false) {
     searchActivatedFlag = true;
     startGoogleCustomSearch();
     }
  }
  else
  {
    toggles[i].style.display = 'none';
    hotspots[i].style.backgroundColor = '#f5f5f5';
    hotspots[i].style.borderBottom = '1px solid #d5d5d5';
    hotspots[i].style.borderRight = '1px solid #f5f5f5';
    hotspots[i].style.borderLeft = '1px solid #f5f5f5';
    document.getElementById('extended_menu_outer').style.display = 'none';
    document.getElementById('search_wrapper_outer').style.style.display = 'none';
  }
}

function hideall()
{
  document.getElementById('extended_menu_outer').style.display = 'none';
  document.getElementById('search_wrapper_outer').style.display = 'none';
  for (var i = 0; i < toggles.length; i++)
  {
    toggles[i].style.display = 'none';
    hotspots[i].style.backgroundColor = '#f5f5f5';
    hotspots[i].style.borderBottom = ' 1px solid #d5d5d5';
    hotspots[i].style.borderRight = ' 1px solid #f5f5f5';
    hotspots[i].style.borderLeft = ' 1px solid #f5f5f5';
  }
}
